package hmm;

import static org.junit.Assert.*;

import org.junit.Test;

public class HmmTest
{
	int M = 20; /*symbols*/
	int N = 50; /*states*/
	
	@Test
	public void bigTest()
	{
		
		int T_MAX = 250;
		
		int obs[] = new int[T_MAX];
		double[] iProb = new double[N];
		double[][] tProb = new double[N][N];
		double[][] eProb = new double[M][N];
		int[] vocab = new int[M]; 
		
		for(int i=0; i < M; ++i)
			vocab[i] = i;
		
		init(obs, iProb, tProb, eProb);
		
		Hmm h = new Hmm(N, vocab, obs, iProb, tProb, eProb); 
		System.out.println("Big Test");
		h.forward();
		System.out.println(h.findFinalProbability()/Math.log(10));
		h.backward();
		System.out.println(h.findFinalProbabilityBackward()/Math.log(10));
		
	}
	
	@Test
	public void smallTest()
	{
		int T_MAX = 40;
		int obs[] = new int[T_MAX];
		double[] iProb = new double[N];
		double[][] tProb = new double[N][N];
		double[][] eProb = new double[M][N];
		int[] vocab = new int[M]; 
		
		for(int i=0; i < M; ++i)
			vocab[i] = i;
		
		init(obs, iProb, tProb, eProb);
		
		Hmm h = new Hmm(N, vocab, obs, iProb, tProb, eProb); 
		System.out.println("Small Test");
		h.forward();
		System.out.println(h.findFinalProbability()/Math.log(10));
		h.backward();
		System.out.println(h.findFinalProbabilityBackward()/Math.log(10));
	}
	
	public void init(int[] observations, double[] initProb, double[][] transitProb, double[][] emitProb)
	{
		double p = 1.0 / N;
		for(int i=0; i < N; ++i)
		{
			initProb[i] = p; 
			for(int j=0; j < N; ++j)
				transitProb[i][j] = p; 
		}
		
		p = 1.0 / M; 
		for(int i=0; i < M; ++i)
		{	
			observations[i] = i;  
			for(int j=0; j < N; ++j)
				emitProb[i][j] = p;
		}
	}

}
