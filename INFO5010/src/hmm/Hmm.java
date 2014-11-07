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
	
	private double alpha[][];	//At time [t], the probability of being in state [i] and observing X(1)..X(t)
	private double beta[][];	//At time [t], given that the state is [i], the probability of observing X(t+1)..X(T)
	private double gamma[][];	//At time [t], the probability of being in state [i]
	private int observations[];
	private	int totalTimeSteps;
	
	private double scaleFactor[];
	
	public Hmm(int numStates, int[] vocab, int[] observations, double[] initProb, double[][] transitProb, double[][] emissionProb)
	{
		N = numStates;
		states = new int[numStates];
		this.initProb = initProb;
		//Better to abstract this out
		this.transitProb = transitProb; //From state in the column to state in the row i.e. transitProb[to][from]
		this.emissionProb = emissionProb;
		this.observations = observations;
		totalTimeSteps = observations.length;
		
		alpha = new double[totalTimeSteps][numStates];
		beta = new double[totalTimeSteps][numStates];
		gamma = new double[totalTimeSteps][numStates];
		scaleFactor = new double[totalTimeSteps];
	}
	
	
	public void checkModelValidity() throws Exception
	{
		double err = 1e-5;
		double total = 0;
		
		//Number of States has to be more than 0
		if(N <= 0)
			throw new Exception("Number of states must be greater than 0");
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
		
		//
	}
	
	public void forward()
	{
		//Base case, t = 0.  
		scaleFactor[0] = 0.0;
		for(int i = 0; i < N; ++i)
		{
			alpha[0][i] = emissionProb[observations[0]][i] * initProb[i];
			scaleFactor[0] += alpha[0][i];
		}
		//Apply scaling
		scale(alpha[0], scaleFactor[0]);
		//Calculate the remaining time steps
		for(int t = 1; t < totalTimeSteps; ++t)
		{
			scaleFactor[t] = 0.0;
			for(int i = 0; i < N; ++i)
			{
				double sum = 0.0;
				for(int j = 0; j < N; ++j)
					sum += transitProb[i][j] * alpha[t-1][j];
				alpha[t][i] = emissionProb[observations[t]][i] * sum;
				scaleFactor[t] += alpha[t][i];				
			}
			//Apply Scaling
			scale(alpha[t], scaleFactor[t]);
		}		
	}
	
	public double findFinalProbability()
	{
		//Return the probability that the given HMM generated the given sequence
		double result = 0.0;
		for(int i = 0; i < N; ++i)
			result += alpha[totalTimeSteps - 1][i];
		//Reverse the scaling NB; log e
		result = Math.log(result);
		for(int i = 0; i < totalTimeSteps; ++i)
			result += Math.log(scaleFactor[i]);
		return result;
	}
	
	//???
	public double findFinalProbabilityBackward()
	{
		//Return the probability that the given HMM generated the given sequence
		double result = 0.0;
		for(int i = 0; i < N; ++i)
			result += beta[0][i];
		//Reverse the scaling NB; log e
		result = Math.log(result);
		for(int i = 1; i < totalTimeSteps; ++i)
			result += Math.log(scaleFactor[i]);
		return result;
	}

	//APPLY SCALING? using constants from alpha is suggested ? ERROROR
	public void backward()
	{
		//Base Case, t = T. Probability of emitting an empty sequence
		for(int i=0; i < N; ++i)
			beta[totalTimeSteps - 1][i] = 1.0;
		//Calculate probability for previous time-steps
		for(int t = totalTimeSteps - 2; t >= 0; --t)
		{
			for(int from = 0; from < N; ++from)
			{
				double sum = 0.0;
				for(int to = 0; to < N; ++to)
				{
					sum += transitProb[to][from] * emissionProb[observations[t+1]][to] * beta[t+1][to];
				}
				beta[t][from] = sum;
				//Apply scaling NB: one step ahead
				scale(beta[t], scaleFactor[t+1]);
			}
		}
	}
	
	public void computeGamma()
	{
		
		for(int t=0; t < totalTimeSteps; ++t)
		{
			double sum = 0.0;
			for(int i=0; i < N; ++i)
			{
				gamma[t][i] = alpha[t][i] * beta[t][i];
				sum += gamma[t][i];
			}
			for(int i=0; i < N; ++i)
				gamma[t][i] /= sum;
		}
	}
	
	public int[] localStateEstimates()
	{
		int[] optimalStates = new int[totalTimeSteps];
		
		for(int t=0; t < totalTimeSteps; ++t)
		{
			double bestProb = gamma[t][0];
			optimalStates[t] = 0;
			for(int i=0; i < N; ++i)
			{
				if(bestProb < gamma[t][i])
				{
					bestProb = gamma[t][i];
					optimalStates[t] = i;
				}
			}
		}
		return optimalStates;
	}
	
	
	public int[] viterbi() throws Exception
	{
		double[][] delta = new double[totalTimeSteps][N];
		int[][] psi = new int[totalTimeSteps][N]; 
		
		for(int i=0; i < N; ++i)
		{
			delta[0][i] = emissionProb[observations[0]][i] * initProb[i];
		}
		
		for(int t=1; t < totalTimeSteps ; ++t)
		{
			double scaleFactor = 0;
			for(int to = 0; to < N; ++to)
			{
				int bestState = 0;
				double bestStateProb = 0.0;
				for(int from = 0; from < N; ++from)
				{
					double res = transitProb[to][from] * delta[t-1][from];
					if(bestStateProb < res)
					{
						bestStateProb = res;
						bestState = from;
					}					
				} 
				psi[t][to] = bestState;
				delta[t][to] = emissionProb[observations[t]][to] * bestStateProb; 
				scaleFactor += delta[t][to];
			}
			scale(delta[t], scaleFactor);
		}
		
		int[] optimalStates = new int[totalTimeSteps];
		//Find the Final optimal State
		optimalStates[totalTimeSteps - 1] = argMaxState(delta[totalTimeSteps - 1]);
		//Backtrack to find preceding states
		for(int t = totalTimeSteps - 1; t > 0; --t)
			optimalStates[t-1] = psi[t][optimalStates[t]];
		
		return optimalStates;
	}
	
	public int argMaxState(double x[]) throws Exception
	{
		int best;
		double bestVal; 
		
		if(x.length == 0)
			throw new Exception("argMax requires a non-empty array");
		
		best = 0;
		bestVal = x[0];
		for(int i = 1; i < N; ++i)
		{
			if(bestVal < x[i])
			{
				best = i;
				bestVal = x[i];
			}
		}
		return best;
	}
	
	public int argMin(double x[]) throws Exception
	{
		int best;
		double bestVal; 
		
		if(x.length == 0)
			throw new Exception("argMax requires a non-empty array");
		
		best = 0;
		bestVal = x[0];
		for(int i = 1; i < N; ++i)
		{
			if(bestVal > x[i])
			{
				best = i;
				bestVal = x[i];
			}
		}
		return best;
	}
	
	public double dotProduct(double[] x, double[] y) throws Exception
	{
		double s = 0.0;
		if(x.length != y.length)
			throw new Exception("Vector size must be equal for dot product.");
		for(int i=0; i < x.length; ++i)
			s += x[i] * y[i];
		return s;
	}
	
	public double sum(double[] x)
	{
		double s = 0.0;
		for(int i=0; i < x.length; ++i)
			s += x[i];
		return s;
	}
	
	public void scale(double x[], double c)
	{
		for(int i=0; i < x.length; ++i)
			x[i] /= c;
	}
	
	public double[] elemProduct(double[] x, double[] y) throws Exception
	{
		double[] res =  new double[x.length];
		if(x.length != y.length)
			throw new Exception("Vector size must be equal for element product.");
		for(int i=0; i < x.length; ++i)
			res[i] += x[i] * y[i];
		return res;
	}
	


}
		
