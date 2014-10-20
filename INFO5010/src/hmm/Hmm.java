package hmm;

/** 
 *  Implements an Hidden Markov Model.
 *  	Contains : 
 *  		states: 
 *  				Each state has a given probability of being initial and an associated observation probability function. 
 *  				Each state is indexed from 0 to n-1, where n is the total number of states. 
 *   		transition probabilities: 
 */


public class Hmm 
{
	private int N;
	private int states[]; 
	private double initProb[]; 
	private double transitProb[][];
	private double emissionProb[][]; 
	
	private double alpha[][];
	private double beta[][];
	private int observations[];
	private	int totalTimeSteps;
	
	public Hmm(int numStates, String[] vocab, int[] observations)
	{
		N = numStates;
		states = new int[numStates];
		initProb = new double[numStates];
		transitProb = new double[numStates][numStates];
		emissionProb = new double[vocab.length][numStates];
		this.observations = observations;
		totalTimeSteps = observations.length;
		
		alpha = new double[totalTimeSteps][numStates];
		beta = new double[totalTimeSteps][numStates];
		
	}
	
	
	public void checkModelValidity() throws Exception
	{
		double err = 1e-5;
		double total = 0;
		
		//Sum of all initial Probabilities = 1
		for(int i=0; i < N; ++i)
			total += initProb[i];
		if(!(total < 1.0 + err && total > 1.0 - err))
			throw new Exception("Initital Probabilities must sum to 1.0");
		
		//State transition probabilities must sum to 1 for any state
		for(int j=0; j < N; ++j)
		{
			total = 0;
			for(int i=0; i < N; ++i)
				total += transitProb[i][j];
			if(!(total < 1.0 + err && total > 1.0 - err))
				throw new Exception("Transition Probabilities for each state must sum to 1.0");
		}
	}
	
	public double forward()
	{
		//Base case, t = 0.  
		for(int i = 0; i < N; ++i)
			alpha[0][i] = emissionProb[observations[0]][i] * initProb[i];
		
		//Calculate the remaining time steps
		for(int t = 1; t < totalTimeSteps; ++t)
		{
			for(int i = 0; i < N; ++i)
			{
				double sum = 0.0;
				for(int j = 0; j < N; ++j)
					sum += transitProb[i][j] * alpha[t-1][j];
				alpha[t][i] = emissionProb[observations[t]][i] * sum;
			}
		}
		
		//Return the probability that the given HMM generated the given sequence
		double result = 0.0;
		for(int i = 0; i < N; ++i)
			result += alpha[totalTimeSteps - 1][i];
		return result; 
		
	}


	public double backward()
	{
		//Base Case, t = T. Probability of emitting an empty sequence
		for(int i=0; i < N; ++i)
			beta[totalTimeSteps - 1][i] = 1.0;
		
		//Calculate probability for previous time-steps
		for(int t = totalTimeSteps - 1; t > 0; --t)
		{
			for(int i = 0; i < N; ++i)
			{
				double sum = 0.0;
				for(int j = 0; j < N; ++j)
				{
					sum += emissionProb[observations[t][i] * transitProb[]
				}
		}
		
	}
	}
}
