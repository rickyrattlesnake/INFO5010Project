package sequence;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Set;

/**
 * Represents a sequence of data, given an alphabet
 * @author Ricky
 *
 */ 
public class Sequence
{
	private Alphabet alphabet;
	private ArrayList<Integer> sequence;
  
	public Sequence(Alphabet alphabet, String sequenceString, String delimiter)
	{
		this.alphabet = alphabet;
		sequence = new ArrayList<Integer>();
		
		/*Check for empty Sequence*/
		if(sequenceString.length() > 0)
		{
			for(String sym : sequenceString.split(delimiter))
			{
				sequence.add(alphabet.getInt(sym));
			}
		}
	}
	
	public Sequence(Alphabet alphabet, String sequenceString)
	{
		this(alphabet, sequenceString, "");
	}
	
	
	/**
	 * Inserts the given motif into a random position in the sequence 
	 * @param motif
	 * @return start index of the motif in the sequence
	 */
	public int insertMotif(Sequence motif)
	{
		try{
			if(!this.alphabet.equals(motif.alphabet))
				throw new Exception("Alphabet of inserted motif must be the same as the Sequence");
		}catch(Exception e){
			e.printStackTrace();
			System.exit(1);
		}
		
		Random gen = new Random();
		int insertionIndex = gen.nextInt(getSize());
		sequence.addAll(insertionIndex, motif.sequence);
		
		return insertionIndex;
	}
  
	
	/**
	 * Inserts a mutation given the probability of mutation for each position in the sequence
	 * mutation probability is given by the underlying distribution of the alphabet
	 * @param probMutation
	 */
	public void mutate(double probMutation)
	{
		Random gen = new Random();
		for(int i=0; i < getSize(); ++i)
		{
			double randomVar = gen.nextDouble();
			if(randomVar < probMutation)
			{
				sequence.set(i, alphabet.getRandomInt());
			}
		}
	}
	
	public static Sequence generateRandomSequence(Alphabet alphabet, int size)
	{
		Sequence output = new Sequence(alphabet, "");
		
		for(int i=0; i < size; ++i)
		{
			output.append(alphabet.getRandomSymbol());
		}
		
		return output;
	}
	
	/**
	 * Returns a subsequence from [start, end)
	 * @param start
	 * @param end
	 * @return
	 */
	public Sequence getSubsequence(int start, int end)
	{
		//TODO: need to make this more efficient
		Sequence subSeq = new Sequence(alphabet, "");
		for(int i=start; i < end; ++i)
		{
			subSeq.append(this.getPosition(i));
		}
		return subSeq;
	}
	
	public void append(String sym)
	{
		try{
			if(!alphabet.contains(sym))
				throw new Exception("Symbol to append is not in the Alphabet");
		}catch(Exception e){
			e.printStackTrace();
			System.exit(1);
		}
		sequence.add(alphabet.getInt(sym));
	}
	
	public List<Integer> sequenceToInt()
	{
		return sequence;
	}
	
	public Integer getIntAtPosition(int i)
	{
		return sequence.get(i);
	}
	
	public String getPosition(int i)
	{
		return alphabet.getSymbol(getIntAtPosition(i));
	}
	
	
	public int getSize()
	{
		return sequence.size();
	}
	
	public Alphabet getAlphabet()
	{
		return alphabet;
	}
	
	@Override
	public String toString()
	{
		StringBuilder output = new StringBuilder();
		
		for(Integer i : sequence)
		{
			output.append(alphabet.getSymbol(i));
		}
		
		return output.toString();
	}

	/**
	 * Returns the projection of a template onto the l-mer starting at the given position.
	 * @param start
	 * @param template
	 * @return String projected symbols in order of index
	 */
	public String getProjection(int start, Set<Integer> template)
	{
		List<Integer> sortedTemplate = new ArrayList<Integer>(template);
		Collections.sort(sortedTemplate);
		
		String projection = "";
		for(Integer pos : sortedTemplate)
		{
			projection += getPosition(pos + start);
		}
		
		return projection;
		
	}
	
	/**
	 * Creates a copy of the current sequence
	 * @return
	 */
	public Sequence copy()
	{
		Sequence copy = new Sequence(alphabet, this.toString());
		return copy;
	}
}
