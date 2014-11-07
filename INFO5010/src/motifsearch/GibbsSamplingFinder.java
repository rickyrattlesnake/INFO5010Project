package motifsearch;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import scoring.ExpectationScore;
import scoring.RelativeInformationScore;
import scoring.Score;
import sequence.Alphabet;
import sequence.Profile;
import sequence.Sequence;

public class GibbsSamplingFinder extends Finder
{
	private double optimizationThreshold; //The convergence test threshold
	
	public GibbsSamplingFinder(Alphabet alphabet, List<Sequence> seqList,
			int motifLength, double optimizationThreshold, Score scorer)
	{
		super(alphabet, seqList, motifLength, scorer);
		this.optimizationThreshold = optimizationThreshold;
	}
	
	public GibbsSamplingFinder(Alphabet alphabet, List<Sequence> seqList,
			int motifLength, double optimizationThreshold)
	{
		super(alphabet, seqList, motifLength, new RelativeInformationScore());
		this.optimizationThreshold = optimizationThreshold;
	}
	
	/**
	 * Convergence is determined by TODO
	 * @param alphabet
	 * @param seqList
	 * @param motifLength
	 * @return
	 */
	@Override
	public Sequence findMotifs()
	{
		printAlgorithmStart("Gibbs Sampling Finder");
		Random gen = new Random();
		boolean scoreChanged = true;
		
		//Randomly select starting positions
		currentProfile.generateRandomAlignment();
		
		double currentProfileScore = scorer.calculateScore(currentProfile);
		int iters = 0;
		while(scoreChanged)
		{				
			scoreChanged = false;
			ArrayList<Sequence> unoptimizedSequences = new ArrayList<Sequence>(seqList);
			
			/* Iterate randomly through the sequence list and optimize
			 * each alignment whilst checking the profile score for each optimization. 
			 */
			while(!unoptimizedSequences.isEmpty())
			{
				if(iters % 100 == 0)
					printIterationInfo(iters, currentProfileScore, currentProfile.alignmentsToString());
				
				//Randomly select one sequence from the unoptimized sequences
				int selection = gen.nextInt(unoptimizedSequences.size());
				Sequence selectedSeq = unoptimizedSequences.get(selection);
				//Create a new profile for sequences without the selected sequence
				List<Sequence> culledSeqList = new ArrayList<Sequence>(seqList);
				culledSeqList.remove(selectedSeq);
				Map<Sequence, Integer> culledAlignments = currentProfile.getAlignmentStarts();				
				culledAlignments.remove(selectedSeq);	
				Profile culledProfile = new Profile(alphabet, culledSeqList, motifLength, culledAlignments);
				
				//For each position for the selected sequence find the distribution of scores
				double[] scores = culledProfile.scoreAllLmers(selectedSeq, scorer);
				
				/* Pick a random new alignment start index based on score distribution. 
				 * Since scores can be negative, we shift all scores by the lowest score,
				 * to maintain distribution properties and using cumulative distributions we choose 
				 * the new index
				 */
				double minScore = Double.MAX_VALUE;
				for(int i=0; i < scores.length; ++i)
				{
					minScore = (scores[i] < minScore) ? scores[i] : minScore;
				}
				
				double cumulativeScore = 0.0;
				for(int i=0; i < scores.length; ++i)
				{
					scores[i] += Math.abs(minScore);
					cumulativeScore += scores[i];
				}
				
				double threshold = gen.nextDouble() * cumulativeScore;
				
				int newStart = 0;
				for(int i=0; i < scores.length; ++i)
				{
					if(threshold < scores[i])
					{
						newStart = i;
						break;
					}
					else
					{
						threshold -= scores[i];
					}
				}
				
				//Update the new alignment for selected sequence
				double oldScore = scorer.calculateScore(currentProfile);
				int oldStart = currentProfile.getAlignmentStarts().get(selectedSeq);
				currentProfile.updateAlignmentStart(selectedSeq, newStart); 
				double newScore = scorer.calculateScore(currentProfile);
				
				//Check for changes 
				if((newScore - oldScore) > optimizationThreshold)
				{
					//If profile score has improved, restart random optimization
					scoreChanged = true;
					break;
				}
				else
				{
					//If profile score has not improved, select another alignment to optimize 
					//First reverse the change made by the start change
					currentProfile.updateAlignmentStart(selectedSeq, oldStart); 
					unoptimizedSequences.remove(selection);					
				}
			}
			currentProfileScore = scorer.calculateScore(currentProfile);
			iters++;
		}
		printAlgorithmEnd(currentProfileScore, currentProfile.getConsensus());
		return currentProfile.getConsensus();
	}
}
