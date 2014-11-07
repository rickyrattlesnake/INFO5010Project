package scoring;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import sequence.Alphabet;
import sequence.Profile;
import sequence.Sequence;

public class ScoreTest
{

	@Test
	public void test()
	{
		Alphabet alph = new Alphabet("ACGT", "");
		Sequence s1 = new Sequence(alph, "CAGGACACTTTCTAAACTGCCTAACTAAGATGCGTGCCCTTCGATTTTCAGGCTGTT");
		Sequence s2 = new Sequence(alph, "GGAGGATACTATCAGTATTATACACCAGCGCTTCTTTCGGATTTTTAGAGCCCTGTG");
		Sequence s3 = new Sequence(alph, "TAAGTAGCGTGGTCAACAACGTTGGCTAACAGGAAGGGCCAAAATTATTAGTGGAAG");
		Sequence s4 = new Sequence(alph, "GAGAAACGGACATGGTGTACATTGGTCGGCTGTGGAATTGTATGCTCAGGTCTGGCT");
		Sequence s5 = new Sequence(alph, "GTGAAATTCAACTCAGGTATGATTGACGTCCGGCCTGTGGCGGAGCGGCTCAGGGCA");

		List<Sequence> slist = new ArrayList<Sequence>();
		slist.add(s1); 
		slist.add(s2);
		slist.add(s3);
		slist.add(s4);
		slist.add(s5);
		Profile prof = new Profile(alph, slist, 6);
		prof.update();
		
		
		Sequence lmer = new Sequence(alph, "AAACGA");
		Score pscore = new FrequencyScore();
		
		double[] expected = {1e-5, 0.6, 0.4, 1e-5, 0.2, 0.8};
		double exp = 0;
		for(int i=0; i < expected.length; ++i)
		{
			exp += Math.log(expected[i]) / Math.log(2);
		}
		double res = pscore.calculateScore(prof, lmer);

		System.out.println(prof);
		System.out.println(res);
		System.out.println(exp);
		prof.updateAlignmentStart(s1, 4);
		System.out.println(prof);
		System.out.println(prof.getConsensus());
		prof.printPwm();
		System.out.println(pscore.calculateScore(prof));
		
		System.out.println(prof.alignmentsToString());
	}

}
