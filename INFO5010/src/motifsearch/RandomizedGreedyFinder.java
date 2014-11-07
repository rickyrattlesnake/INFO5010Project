package motifsearch;

import java.util.List;
import java.util.Map;

import scoring.RelativeInformationScore;
import scoring.Score;
import sequence.Alphabet;
import sequence.Sequence;

public class RandomizedGreedyFinder extends Finder
{
	private long maxIterations = 50000; 
	private boolean updateEachStep;
	
	/**
	 * Setup for randomized greedy algorithm
	 * @param alphabet
	 * @param seqList
	 * @param motifLength
	 * @param updateEachStep	profile matrix recalculation after all re-alignments or after each re-alignment
	 * @param scorer
	 */
	public RandomizedGreedyFinder(Alphabet alphabet, List<Sequence> seqList,
			int motifLength, boolean updateEachStep, Score scorer)
	{
		super(alphabet, seqList, motifLength, scorer);
		this.updateEachStep = updateEachStep;
	}
	
	/**
	 * Use the default InformationContentScore as a metric
	 * @param alphabet
	 * @param seqList
	 * @param motifLength
	 * @param maxIterations
	 */
	public RandomizedGreedyFinder(Alphabet alphabet, List<Sequence> seqList,
			int motifLength, boolean updateEachStep)
	{
		super(alphabet, seqList, motifLength, new RelativeInformationScore());
		this.updateEachStep = updateEachStep;
	}

	/**
	 * Starts with random motif positions and constructs the profile.
	 * Use the profile to find P-most-probable-lmers in each sequence, and
	 * use these as the new motif positions to construct the new profile. Iterate until 
	 * convergence or max-iterations is reached
	 * @param profile
	 */
	public Sequence findMotifs()
	{		
		printAlgorithmStart("Randomized Greedy Finder");
		
		//Randomly set the motif start positions and update
		currentProfile.generateRandomAlignment();
		
		//optimize the profile by iteratively finding
		//the best l-mers and recalculating profiles
		double bestProfileScore = Double.MIN_VALUE;
		double currentProfileScore = scorer.calculateScore(currentProfile);
		int iters = 0;
		
		while(currentProfileScore > bestProfileScore 
					&& iters < maxIterations)
		{ 
			//printIterationInfo(iters, currentProfileScore, currentProfile.alignmentsToString());
			bestProfileScore = currentProfileScore; 
			
			Map<Sequence, Integer> alignments = currentProfile.getAlignmentStarts();
			for(Sequence seq : seqList)
			{
				//Find the best l-mer
				double[] scores = currentProfile.scoreAllLmers(seq, scorer);
				int bestLmerStart = findMaxIndex(scores);
				//Move the start position to the best l-mer and update if 
				//flag is set by user
				if(updateEachStep)
				{
					currentProfile.updateAlignmentStart(seq, bestLmerStart);
				}
				else
				{
					alignments.put(seq, bestLmerStart);
				}
			}
			//Update after all steps if flag is set by user
			if(!updateEachStep)
			{
				for(Sequence seq : alignments.keySet())
				{
					currentProfile.updateAlignmentStart(seq, alignments.get(seq));
				}
			}
			
			currentProfileScore = scorer.calculateScore(currentProfile);
			iters++;
		}
		
		printAlgorithmEnd(currentProfileScore, currentProfile.getConsensus());
		//The motif is the consensus of the optimized profile
		return currentProfile.getConsensus();		
	}
	
	public int findMaxIndex(double[] x)
	{
		try{
			if(x.length == 0)
				throw new Exception("Empty array error");
		}catch(Exception e){
			e.printStackTrace();
			System.exit(1);
		}
		
		double max = x[0];
		int idx = 0;
		for(int i=1; i < x.length; ++i)
		{
			if(max < x[i])
			{
				max = x[i];
				idx = i;
			}
		}		
		return idx;
	}
}
