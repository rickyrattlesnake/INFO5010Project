package scoring;

import java.util.Arrays;

import sequence.Profile;
import sequence.Sequence;

public class RelativeInformationScore extends Score
{

	public RelativeInformationScore(double pseudoZero)
	{
		super(pseudoZero);
	}
	
	public RelativeInformationScore()
	{
		super(1e-5);
	}
	/**
	 * Calculates the Kullback-Leibler measure for the whole profile Matrix. 
	 * NB: this score will always be >= 0, since the Sum(p * log_2(p/b) of all elements, 
	 * where p = probability of a symbol at a position and b = background probability for the symbol
	 * always provide more or the same level of information as the background.
	 */
	@Override
	public double calculateScore(Profile profile)
	{
		double result = 0;
		for(String symbol : profile.getAlphabet().getAllSymbols())
		{
			for(int pos=0; pos < profile.length(); ++pos)
			{
				result += profile.getPpm(symbol, pos) * profile.getPwm(symbol, pos);
			}
		}
		return result;	
	}
	
	/**
	 * Calculates the Kullback-Leibler measure for the the given l-mer. 
	 * i.e. Sum(p * log_2(p/b)), where p = probability of a symbol at a position, 
	 * and b = background probability for the symbol. This score can be negative as 
	 * we are not summing over the whole matrix.
	 */
	@Override
	public double calculateScore(Profile profile, Sequence lmer)
	{		
		double result = 0;
		for(int pos=0; pos < lmer.getSize(); ++pos)
		{
			result += profile.getPpm(lmer.getPosition(pos), pos) * profile.getPwm(lmer.getPosition(pos), pos);
		}
		return result;
	}
	
	

}
