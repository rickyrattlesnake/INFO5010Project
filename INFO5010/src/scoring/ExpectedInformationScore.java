package scoring;

import sequence.Profile;
import sequence.Sequence;

public class ExpectedInformationScore extends Score
{
	public ExpectedInformationScore(double pseudoZero)
	{
		super(pseudoZero);
	}
	
	public ExpectedInformationScore()
	{
		super(1e-5);
	}
	
	/**
	 * Calculates the Shannon's Information Score for the whole profile matrix. Sum(p * Log_2(p))
	 * @return total score for profile
	 */
	@Override
	public double calculateScore(Profile profile)
	{
		double result = 0;
		for(String symbol : profile.getAlphabet().getAllSymbols())
		{
			for(int pos=0; pos < profile.length(); ++pos)
			{
				double prob = profile.getPpm(symbol, pos);
				
				if(prob < pseudoZero)
				{
					prob = pseudoZero;
				}
				result += prob * (Math.log(prob) / Math.log(2));
			}
		}
		return result;
	}
	
	/**
	 * For the sequence in question calculate the Shannon's Information Score using the given profile
	 */
	public double calculateScore(Profile profile, Sequence lmer)
	{
		double result = 0;
		for(int pos=0; pos < lmer.getSize(); ++pos)
		{
			String symbol = lmer.getPosition(pos);
			double prob = profile.getPpm(symbol, pos);
				
			if(prob < pseudoZero)
			{
				prob = pseudoZero;
			}
			result += prob * (Math.log(prob) / Math.log(2));
		}
		return result;
	}

}
