package motifsearch;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import scoring.Score;
import sequence.Alphabet;
import sequence.Profile;
import sequence.Sequence;

public abstract class Finder
{
	protected Alphabet alphabet;
	protected List<Sequence> seqList; 
	protected int motifLength; 
	protected Profile currentProfile;
	protected Score scorer;
	
	public Finder(Alphabet alphabet, List<Sequence> seqList, int motifLength, Score scorer)
	{
		this.alphabet = alphabet;
		this.seqList = seqList;
		this.motifLength = motifLength;
		currentProfile = new Profile(alphabet, seqList, motifLength);
		this.scorer = scorer;
	}
	
	/**
	 * Runs the finder over multiple trials and returns the motifs corresponding to
	 * the highest scoring trial
	 * @param trials
	 * @return
	 */
	public Sequence runMultiple(int trials)
	{
		Map<Sequence, Double> results = new HashMap<Sequence, Double>();
		
		System.out.println("+++ Starting Multiple Trials +++");
		for(int i=0; i < trials; ++i)
		{
			Sequence motif = findMotifs();
			double score = scorer.calculateScore(currentProfile);
			results.put(motif, score);
			System.out.println(String.format("Trial : %d - Score : %.5f - Motif : %s", i, score, motif.toString()));
		}

		Map.Entry<Sequence, Double> maxEntry = null;

		for (Map.Entry<Sequence, Double> entry : results.entrySet())
		{
		    if (maxEntry == null || entry.getValue().compareTo(maxEntry.getValue()) > 0)
		    {
		        maxEntry = entry;
		    }
		}
		System.out.println(
				String.format("+++ Multiple Trials finished - Max Score : %.5f - Best Motif : %s +++", 
											maxEntry.getValue(), maxEntry.getKey().toString()));
		return maxEntry.getKey();
	}
		
	public abstract Sequence findMotifs();
	
	public Profile getCurrentProfile()
	{
		return currentProfile;
	}
	
	public void printAlgorithmStart(String name)
	{
		System.out.println(String.format("*** Running %s ***", name));
	}
	
	public void printIterationInfo(int iteration, double score, String alignmentArrayString)
	{
		System.out.println(String.format("Iteration : %d - Score : %.5f - Alignments : %s", iteration, score, alignmentArrayString));
	}
	
	public void printAlgorithmEnd(double score, Sequence motif)
	{
		System.out.println(String.format("*** Result : Final Score %.5f - Predicted Motif : %s ***", score, motif.toString()));
	}
	
}
