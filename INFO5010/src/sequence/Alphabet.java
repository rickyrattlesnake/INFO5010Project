package sequence;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

/**
 * Generates and holds an Integer mapping between symbols in the given alphabet.
 * The alphabet cannot be empty, and mapped integers range from [0, size of alphabet) 
 * with no gaps. 
 * @author Ricky
 *
 */
public class Alphabet
{
	private TreeSet<String> symbolSet;				//Allows for multiple character alphabet symbols
	private Map<String, Integer> symbolToIntMap;
	private Map<Integer, String> intToSymbolMap;
	private Map<String, Double> probDistribution; 	//The probability of generating the symbol
	
	/**
	 * Constructs an Alphabet
	 * @param alphabet string array for each symbol
	 * @param probDist background probability of each corresponding symbol
	 */
	public Alphabet(String[] alphabet, double[] probDist)
	{
		this.symbolSet = new TreeSet<String>();
		symbolToIntMap = new HashMap<String, Integer>();
		intToSymbolMap = new HashMap<Integer, String>();
		probDistribution = new HashMap<String, Double>();
		
		/*Check for empty alphabet*/
		if(alphabet.length < 1 || alphabet[0] == "")
		{
			System.err.println("Cannot have an empty Alphabet");
			System.exit(1);
		}
			
		for(int i = 0; i < alphabet.length; ++i)
		{
			String symbol = alphabet[i];
			this.symbolSet.add(symbol);
			//Create a Mapping symbol -> Integer
			symbolToIntMap.put(symbol, i);
			//Create a Mapping Integer -> symbol
			intToSymbolMap.put(i, symbol);
			//Mapping RefInteger -> probability of generating symbol
			probDistribution.put(symbol, probDist[i]);
			
		}
	}
	
	/**
	 * Constructs an Alphabet, assumes a uniform probability distribution
	 * @param alphString
	 * @param delimiter
	 */
	public Alphabet(String alphString, String delimiter)
	{		
		this(alphString.split(delimiter), 
				genUniformDist(alphString.length()));
	}
	
	/**
	 * Constructs an Alphabet with defined probability distribution
	 * @param alphString delimited string of alphabet characters
	 * @param delimiter delimiter character
	 * @param probDist background probability of each corresponding symbol 
	 */
	public Alphabet(String alphString, String delimiter, double[] probDist)
	{
		this(alphString.split(delimiter), probDist);
	}
	
	/**
	 * Returns a random symbol from the alphabet using the underlying probability distribution
	 * @return symbol
	 */
	public String getRandomSymbol()
	{
		Random gen = new Random();
		double randomVar = gen.nextDouble();
		//Test random variable against the cumulative probability distribution
		double cummulativeProb = 0.0;
		for(String sym : symbolSet)
		{
			cummulativeProb += getProbability(sym);
			if(randomVar < cummulativeProb)
			{
				return sym;
			} 
		};
		//Case where inequality check failed due to floating point inaccuracy
		return symbolSet.last();
	}
	
	/**
	 * Helper method used to generate uniform distribution of a given size
	 * @param size
	 * @return probDist
	 */
	private static double[] genUniformDist(int size)
	{
		double[] probDist = new double[size];
		double prob = 1.0 / size;
		for(int i=0; i < size; ++i)
		{
			probDist[i] = prob;
		}
		return probDist;
	}
	
	/**
	 * Constructs a list of all possible words of a given length for 
	 * this alphabet
	 * @param output list to append with words
	 * @param builder used to build the words
	 * @param length length of output word
	 * @param position parameter for recursion (start at 0)
	 */
	public void getAllPossibleSequences(List<String> output, String[] builder, int length, int position)
	{
		if(position == length)
		{
			//Add the string to the output list
			output.add(String.join("", builder));
			return;
		}
		
		for(String symbol : symbolSet)
		{
			builder[position] = symbol; 
			getAllPossibleSequences(output, builder, length, position + 1);
		}
	}
	
	public int getInt(String symbol)
	{
		return symbolToIntMap.get(symbol);		
	}
	
	public String getSymbol(int refInt)
	{
		return intToSymbolMap.get(refInt);
	}
	
	public int getSize()
	{
		return symbolSet.size();
	}
	
	public boolean contains(String symbol)
	{
		return symbolSet.contains(symbol);
	}
	
	public Set<String> getAllSymbols()
	{
		return symbolSet;
	}
	
	public double getProbability(String sym)
	{
		return probDistribution.get(sym);
	}
	
	public int getRandomInt()
	{
		return getInt(getRandomSymbol());
	}
	
	@Override
	public boolean equals(Object other)
	{
		if(other == null) return false; 
		if(other == this) return true; 
		if(!(other instanceof Alphabet)) return false;
		
		//Two Alphabets are the same if they contain the same symbols in any order		
		return symbolSet.equals(((Alphabet)other).symbolSet);
	}

	public Map<String, Double> getProbabilityDistribution()
	{
		return probDistribution;
	}
	
}
