package sequence;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import scoring.Score;

/**
 * Represents the frequency distributions and matrices of a set of sequences
 * and aligned motifs
 * @author Ricky
 *
 */
public class Profile
{
	private Alphabet alphabet;
	private List<Sequence> sequences; 
	private int length;
	private int[][] positionFrequencyMatrix; 		//frequency of each symbol for each position [symbolMappedInt][position]
	private double[][] positionProbabilityMatrix;	//empirical probability of each symbol for each position [symbolMappedInt][position]
	private double[][] positionWeightMatrix; 		//Base 2 log weighted probability of each symbol for each position [symbolMappedInt][position]
	private Map<Sequence, Integer> alignments; 		//Tracks the start of the motif/pattern
	private Map<String, Double> backgroundModel;
	private double DELTA;							//pseudo-zero for avoiding overflow errors when taking logs, and for sampling error
	
	
	public Profile(Alphabet alphabet, List<Sequence> seqList, int length)
	{
		this.alphabet = alphabet;
		//Check if sequences have the same alphabet
		//If no sequence list is given then initialize as an empty list 
		if(seqList == null)
		{
			sequences = new ArrayList<Sequence>();
		}
		else
		{
			for(Sequence s : seqList)
			{
				if(!alphabet.equals(s.getAlphabet()))
				{
					System.err.println("Profiler requires sequences with the same alphabet");
					System.exit(1);
				};
			}
			sequences = seqList;
		}
		
		this.length = length;
		positionFrequencyMatrix = new int[alphabet.getSize()][length];
		positionProbabilityMatrix = new double[alphabet.getSize()][length];
		positionWeightMatrix = new double[alphabet.getSize()][length];
		alignments = new HashMap<Sequence, Integer>();
		
		//Inititate the alignment pointers to zero
		for(Sequence seq : sequences)
		{
			alignments.put(seq, 0);
		}
		
		DELTA = 1.0 / (10 * sequences.size()); 				//Rule of thumb
		
		//Set the background model
		backgroundModel = getEmpiricalBackgroundModel();
		
		//Update the profiles as initialization
		update();		
	}
	
	/**
	 * Constructor with given alignment, makes a copy of the alignment map.
	 * @param alphabet
	 * @param seqList
	 * @param length
	 * @param alignments
	 */
	public Profile(Alphabet alphabet, List<Sequence> seqList, int length, Map<Sequence, Integer> alignments)
	{
		this(alphabet, seqList, length);
		this.alignments = new HashMap<Sequence, Integer>(alignments);
		update();
	}
	
	/**
	 * update the matrix frequencies and probabilities for the current alignment
	 */
	public void update()
	{
		//Clear the profile matrices
		for(int i=0; i < height(); ++i)
		{
			Arrays.fill(positionFrequencyMatrix[i], 0);
			Arrays.fill(positionProbabilityMatrix[i], 0.0);
			Arrays.fill(positionWeightMatrix[i], 0.0);
		}
		
		//Load the frequency data from the sequences
		for(Sequence seq : sequences)
		{
			for(int pos=0; pos < length; ++pos)
			{
				int startMarker = alignments.get(seq);
				int lmerPos = startMarker + pos;
				String symbol = seq.getPosition(lmerPos);
				int previousFreq = positionFrequencyMatrix[alphabet.getInt(symbol)][pos];
				modifyProfileMatrices(symbol, pos, previousFreq + 1);
			}
		}
	}
	
	
	/**
	 * simultaneously update all three profile matrices by changing the frequecny at one position
	 * @param symbol
	 * @param position
	 * @param newFreq
	 */
	private void modifyProfileMatrices(String symbol, int position, int newFreq)
	{
		try{
			if(newFreq < 0)
				throw new Exception("Frequency cannot be less than 0");
		}catch(Exception e){
			e.printStackTrace();
			System.exit(1);
		}
		int mappedInt = alphabet.getInt(symbol);
		positionFrequencyMatrix[mappedInt][position] = newFreq;
		positionProbabilityMatrix[mappedInt][position] = ((double)newFreq) / getSequenceCount();
		double weightedProb = positionProbabilityMatrix[mappedInt][position] / backgroundModel.get(symbol);
		if(weightedProb < DELTA)
		{
			weightedProb = DELTA;
		}
		positionWeightMatrix[mappedInt][position] = Math.log(weightedProb) / Math.log(2);
	}
	
	/**
	 * Calculates the consensus motif for the current profile
	 * @return
	 */
	public Sequence getConsensus()
	{
		Sequence consensus = new Sequence(alphabet, "");
		for(int pos=0; pos < length; ++pos)
		{
			String maxSymbol = "";
			double maxProb = Double.MIN_VALUE;
			for(String symbol : alphabet.getAllSymbols())
			{
				double prob = getPpm(symbol, pos);
				if(prob > maxProb)
				{
					maxSymbol = symbol;
					maxProb = prob;
				}
			}
			consensus.append(maxSymbol);
		}
		
		return consensus;
	}
	
	/**
	 * Changes the start of the alignment for a given sequence in the profile 
	 * and performs an efficient update of the profile matrices
	 * @param seq sequence to update
	 * @param newStartPos new alignment start
	 */
	public void updateAlignmentStart(Sequence seq, int newStartPos)
	{
		try{
			if(!isSequenceInProfile(seq))
				throw new Exception("Sequence not in the profile");
			if(newStartPos < 0 || newStartPos >= seq.getSize() - length + 1)
				throw new Exception("New Start position out of bounds");
		}catch(Exception e){
			e.printStackTrace();
			System.exit(1);
		}
		
		int oldStartPos = alignments.get(seq);
		alignments.put(seq, newStartPos);
				
		for(int i = 0; i < length; ++i)
		{
			//Remove frequency contributions from the old alignment 
			String symbol = seq.getPosition(oldStartPos + i);
			int oldFreq = positionFrequencyMatrix[alphabet.getInt(symbol)][i];
			if(oldFreq > 0)
				modifyProfileMatrices(symbol, i, oldFreq - 1);
			//and add contribution from new alignment
			symbol = seq.getPosition(newStartPos + i);
			oldFreq = positionFrequencyMatrix[alphabet.getInt(symbol)][i];
			modifyProfileMatrices(symbol, i, oldFreq + 1);
		}		
	}
	
	/**
	 * Finds the score for all l-mers from [0, N - length + 1). 
	 * This method is an optimization, to avoid creating subsequences externally 
	 * for finding the highest scoring l-mer where l is the motif length.
	 * @param seq
	 * @param scorer
	 * @return scores of all l-mes [0, n - length + 1)
	 */
	public double[] scoreAllLmers(Sequence seq, Score scorer)
	{
		double[] result = new double[seq.getSize() - length + 1];
		for(int pos=0; pos < result.length; ++pos)
		{
			result[pos] = scorer.calculateScore(this, seq.getSubsequence(pos, pos + length));
		}
		return result;
	}
	
	/**
	 * Generate a random alignment and update the profile matrices
	 */
	public void generateRandomAlignment()
	{
		//Randomly set the motif start positions
		Random gen = new Random();
		for(Sequence s : sequences)
		{
			updateAlignmentStart(s, gen.nextInt(s.getSize() - length + 1));
		}
	}
	
	/**
	 * Find the background probabilities for each symbol
	 * in the profile's sequences
	 * @return
	 */
	public Map<String, Double> getEmpiricalBackgroundModel()
	{
		Map<String, Double> symbolProb = new HashMap<String, Double>();
		
		for(Sequence seq : sequences)
		{
			for(int i=0; i < seq.getSize(); ++i)
			{
				String symbol = seq.getPosition(i);
				if(!symbolProb.containsKey(symbol))
					symbolProb.put(symbol, 0.0);
				symbolProb.put(symbol, symbolProb.get(symbol) + 1);
			}
		}
		
		//Calculate Probabilities
		int total = getTotalSymbolCount();
		for(String sym : symbolProb.keySet())
			symbolProb.put(sym, symbolProb.get(sym) / total);
			
		return symbolProb;
	}
	
	/**
	 * Find the background probabilities for each symbol
	 * in the profile's alphabet (exogenous variable)
	 * @return
	 */	
	public Map<String, Double> getBackgroundModel()
	{
		return alphabet.getProbabilityDistribution();
	}
	
	
	
	/**
	 * Returns the frequency of symbols in the profile at the given symbol and position
	 * @param symbol
	 * @param position
	 * @return
	 */
	public int getPfm(String symbol, int position)
	{
		return positionFrequencyMatrix[alphabet.getInt(symbol)][position];
	}
		
	/**
	 * Returns the probability of the symbol occuring at the specified position
	 * @param symbol
	 * @param position
	 * @return
	 */
	public double getPpm(String symbol, int position)
	{
		return positionProbabilityMatrix[alphabet.getInt(symbol)][position];
	}
	
	/**
	 * Returns the information content of the symbol at the position, 
	 * weighted by the background model
	 * @param symbol
	 * @param position
	 * @return
	 */
	public double getPwm(String symbol, int position)
	{
		return positionWeightMatrix[alphabet.getInt(symbol)][position];
	}
	
	
	/**
	 * Returns a copy of the alignment starts
	 * @return
	 */
	public Map<Sequence, Integer> getAlignmentStarts()
	{
		return new HashMap<Sequence, Integer>(alignments);
	}
	
	/**
	 * Prints the position weight matrix
	 */
	public void printPwm()
	{
		StringBuilder output = new StringBuilder();
		for(int a=0; a < alphabet.getSize(); ++a)
		{
			output.append(String.format("%s | ", alphabet.getSymbol(a)));
			for(int pos=0; pos < length; ++pos)
			{
				output.append(String.format("%.4f ", positionWeightMatrix[a][pos]));
			}
			output.append("\n");			
		}
		System.out.println(output.toString());
	}
	
	/**
	 * Converts the Alignment Map into a string representation of an array
	 * in the order given by the sequence list
	 */
	public String alignmentsToString()
	{
		StringBuilder output = new StringBuilder();
		output.append("[");
		for(Sequence seq : sequences)
		{
			output.append(alignments.get(seq) + ", ");
		}
		output.replace(output.length()-2, output.length(), "");
		output.append("]");
		
		return output.toString();
	}
	
	
	/**
	 * Prints the position probability matrix
	 */
	@Override
	public String toString()
	{
		StringBuilder output = new StringBuilder();
		for(int a=0; a < alphabet.getSize(); ++a)
		{
			output.append(String.format("%s | ", alphabet.getSymbol(a)));
			for(int pos=0; pos < length; ++pos)
			{
				output.append(String.format("%.4f ", positionProbabilityMatrix[a][pos]));
			}
			output.append("\n");			
		}
		return output.toString();
	}
	
	/**
	 * Height of the profile represents the number of symbols mapped to the matrices
	 * @return
	 */
	public int height()
	{
		return alphabet.getSize();
	}
	
	/**
	 * Length of the profile represents the length of the motif
	 * @return
	 */
	public int length()
	{
		return length;
	}
	
	public Alphabet getAlphabet()
	{
		return alphabet;
	}

	public List<Sequence> getSequenceList()
	{
		return sequences;
	}

	
	
	public int getSequenceCount()
	{
		return sequences.size();
	}
	
	public int getTotalSymbolCount()
	{
		int result = 0;
		for(Sequence s : sequences)
		{
			result += s.getSize();
		}
		return result;
		
	}
	public boolean isSequenceInProfile(Sequence seq)
	{
		return sequences.contains(seq);
	}
	
}
