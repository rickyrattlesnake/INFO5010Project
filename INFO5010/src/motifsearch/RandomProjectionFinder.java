package motifsearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import scoring.RelativeInformationScore;
import scoring.Score;
import sequence.Alphabet;
import sequence.Profile;
import sequence.Sequence;

public class RandomProjectionFinder extends Finder
{
	private int projectionSize;  	//represents k in the original algorithm (i.e. number of positions in the template)
	private int binThreshold;		//determines which bins to select in hash table  
	private int numIterations;		//number of iterations of template selection and hashing to perform
	
	public RandomProjectionFinder(Alphabet alphabet, List<Sequence> seqList,
			int motifLength, int projectionSize, int binThreshold, int numIterations, Score scorer)
	{
		super(alphabet, seqList, motifLength, scorer);
		this.projectionSize = projectionSize;
		this.binThreshold = binThreshold;
		this.numIterations = numIterations;
	}
	
	public RandomProjectionFinder(Alphabet alphabet, List<Sequence> seqList,
			int motifLength, int projectionSize, int binThreshold, int numIterations)
	{
		this(alphabet, seqList, motifLength, projectionSize, binThreshold, numIterations, new RelativeInformationScore());
	}
	
	@Override
	public Sequence findMotifs()
	{
		printAlgorithmStart("Random Projection Finder");
		
		//Create a counter of potential motifs start positions for each sequence
		HashMap<Sequence, int[]> potentialMotifCount = new HashMap<Sequence, int[]>();
		for(Sequence seq : seqList)
		{
			potentialMotifCount.put(seq, new int[seq.getSize()]);		//Sequences sizes may differ
			Arrays.fill(potentialMotifCount.get(seq), 0);				
		}
				
		//Create a hash-table of projection counts
		Map<String, Integer> projectionBins = createProjectionBins(alphabet, projectionSize);
		
		//aggregate data obtained for specified number of iterations
		for(int it=0; it < numIterations; ++it)
		{
			//Reset bin values to zero
			for(String s : projectionBins.keySet())
			{
				projectionBins.put(s, 0);
			}
			
			//Select a random k-l template, where k = projectionSize and l=motifLength
			Set<Integer> template = generateRandomTemplate(projectionSize, motifLength);
			for(Sequence seq : seqList)
			{
				for(int i=0; i < seq.getSize() - motifLength + 1; ++i)
				{
					String projection = seq.getProjection(i, template);
					//Increment the count for the bins
					projectionBins.put(projection, projectionBins.get(projection) + 1);
				}
			}
			
			//Record the motif starts that have a higher than threshold projection count
			for(Sequence seq : seqList)
			{
				for(int i=0; i < seq.getSize() - motifLength + 1; ++i)
				{
					String projection = seq.getProjection(i, template);
					if(projectionBins.get(projection) > binThreshold)
					{
						potentialMotifCount.get(seq)[i] += 1;
					}
				}
			}
		}
		
		//The alignment vector is determined from the highest count position
		//for each sequence in potentialMotifStarts
		//TODO: better method necessary , returns the first max, maybe use expectation scores 
		for(Sequence seq : seqList)
		{
			currentProfile.updateAlignmentStart(seq, findMaxIndex(potentialMotifCount.get(seq)));
		}
		
		double currentProfileScore = scorer.calculateScore(currentProfile);
		printAlgorithmEnd(currentProfileScore, currentProfile.getConsensus());
		return currentProfile.getConsensus();
	}
	
	/**
	 * Creates a map of all projections to counts
	 * @param alphabet
	 * @param keySize
	 * @return
	 */
	private Map<String, Integer> createProjectionBins(Alphabet alphabet, int keySize)
	{
		Map<String, Integer> bins = new HashMap<String, Integer>();
		
		//Get all possible combinations of keySize-word strings 
		List<String> allProjections = new ArrayList<String>();
		String[] builder = new String[keySize];
		alphabet.getAllPossibleSequences(allProjections, builder, keySize, 0);
		
		//Load all possible k-template strings
		for(String projection : allProjections)
		{
			bins.put(projection, 0);
		}
		
		return bins;
	}
	
	/**
	 * Generates a random template from a uniform distribution. 
	 * A template is defined as an set of k positions in the range [0, motifLength)
	 * @param k size of the projection
	 * @param l length of the motif
	 * @return
	 */
	private Set<Integer> generateRandomTemplate(int k, int l)
	{
		Set<Integer> template = new HashSet<Integer>();
		
		//Generate k distinct random integers
		Random gen = new Random();
		while(template.size() != k)
		{
			template.add(gen.nextInt(l));
		}
				
		return template;		
	}
	
	public int findMaxIndex(int x[])
	{
		int maxIndex;
		int maxVal; 
		
		maxIndex = -1;
		maxVal = Integer.MIN_VALUE;
		for(int i = 0; i < x.length; ++i)
		{
			if(maxVal < x[i])
			{
				maxIndex = i;
				maxVal = x[i];
			}
		}
		return maxIndex;
	}
}
