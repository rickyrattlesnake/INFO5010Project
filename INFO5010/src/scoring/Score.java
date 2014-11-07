package scoring;

import sequence.Profile;
import sequence.Sequence;

public abstract class Score
{
	protected double pseudoZero; 
	
	public Score(double pseudoZero)
	{
		this.pseudoZero = pseudoZero;
	}

	/**
	 * Calculate the score for the whole profile, some metrics don't 
	 * implement the method.
	 * @param profile
	 * @return
	 */
	public double calculateScore(Profile profile)
	{
		return 0.0;
	}
	
	/**
	 * Calculate the score for a given l-mer using the given profile
	 * @param profile
	 * @param lmer
	 * @return
	 */
	public abstract double calculateScore(Profile profile, Sequence lmer);
}
