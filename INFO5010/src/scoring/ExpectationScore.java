package scoring;

import sequence.Profile;
import sequence.Sequence;

public class ExpectationScore extends Score
{
	
	public ExpectationScore(double pseudoZero)
	{
		super(pseudoZero);
	}
	
	/**
	 * Provides a default pseudo-zero value
	 */
	public ExpectationScore()
	{
		super(1e-5);
	}
	
	/**
	 * Calculates the sum of the log(probability) of the l-mer appearing, given the profile.
	 * Length of the l-mer is assumed to be the length of the profile. NB: this is a negative score 
	 * but maximization rule still applies.
	 */
	@Override
	public double calculateScore(Profile profile, Sequence lmer)
	{
		double result = 0;
		for(int pos=0; pos < lmer.getSize(); ++pos)
		{
			System.out.println(lmer.getPosition(pos));
			double prob = profile.getPpm(lmer.getPosition(pos), pos);
			if(prob < pseudoZero)
			{
				prob = pseudoZero;
			}
			result += Math.log(prob) / Math.log(2);
		}
		return result;
	}

}
