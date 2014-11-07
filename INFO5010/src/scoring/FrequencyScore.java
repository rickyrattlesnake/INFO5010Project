package scoring;

import sequence.Profile;
import sequence.Sequence;


public class FrequencyScore extends Score
{
	public FrequencyScore(double pseudoZero)
	{
		super(pseudoZero);
	}
	
	public FrequencyScore()
	{
		super(0);
	}

	/**
	 * Calculates a simple score using a sum(frequency) of bases for 
	 * each position of the consensus profile
	 *
	 */
	@Override
	public double calculateScore(Profile profile)
	{
		Sequence consensus = profile.getConsensus();
	
		double result = 0;
		for(String symbol : profile.getAlphabet().getAllSymbols())
		{
			for(int pos=0; pos < profile.length(); ++pos)
			{
				if(symbol == consensus.getPosition(pos))
				{
					result += profile.getPfm(symbol, pos);
				}
			}
		}
		
		return result;
	}

	/**
	 * Calculates the number of same symbols comparing the consensus and the given lmer.
	 * This is usually not a good discriminator since motif length is small and 
	 * the background information at each position is disregarded. 
	 */
	@Override
	public double calculateScore(Profile profile, Sequence lmer)
	{
		Sequence consensus = profile.getConsensus();
		
		double result = 0;
		for(int pos=0; pos < lmer.getSize(); ++pos)
		{
			String symbol = lmer.getPosition(pos);
			if(symbol == consensus.getPosition(pos))
				result += 1; 
		}
		return result;
	}
	
	

}
