package motifsearch;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import scoring.ExpectationScore;
import scoring.ExpectedInformationScore;
import scoring.RelativeInformationScore;
import sequence.Alphabet;
import sequence.Profile;
import sequence.Sequence;

public class FinderTest
{

	@Test
	public void test()
	{
		Alphabet alph = new Alphabet("ACGT", "");
		Sequence s1 = new Sequence(alph, "TATCGGAGTGGCCTGCTCACTTTTTGCCGACAGCAAAAGTATCGTTCTCATTACTGGGCCTATATACCACTTCACTACGAGAATTCGCAGTGAAGCCGGGTACACCA");
		Sequence s2 = new Sequence(alph, "TAGACCTTGCTCACTGACATCGCGGGACCATTGCTCAGAGACTGAATTCGAGGAGTCGTCATTGAGGTGAAACCGTTGTTACAGAGTGTAAGTATGTCCGAAAACAG");
		Sequence s3 = new Sequence(alph, "AAATACCGAAGACTACAACAATCGAGAAGGGCTAGAATTCGCGCGTTATTCAACGTCCTCGGGTAACGAAGTGAGCCCTCCGCCATGTCGACCTGAGCTTAGGCCGC");
		Sequence s4 = new Sequence(alph, "ACTTATCTATTGTGAAGTAGGGACCAAACTCAACATGACCAGTGCGCCCTCTCCACCGGATGAGGAAGGGGCTATCCGAATTAGCAAGAATTCGATATACAAGTATG");
		Sequence s5 = new Sequence(alph, "GAATTCGATCCTTTTTGTGAGTAATCCGATTGTTCCTCCCTCTGCGCAATTTAGGTACTCTCACAGAGTTCGTTTGGCTTATTATAGGTTTGCGTCGAAGATCATTT");

		List<Sequence> slist = new ArrayList<Sequence>();
		slist.add(s1); 
		slist.add(s2);
		slist.add(s3);
		slist.add(s4);
		slist.add(s5);
		Profile expected = new Profile(alph, slist, 7);
		expected.updateAlignmentStart(s1, 80);
		expected.updateAlignmentStart(s2, 43);
		expected.updateAlignmentStart(s3, 34);
		expected.updateAlignmentStart(s4, 87);
		expected.updateAlignmentStart(s5, 0);
		
		System.out.println(expected.getConsensus());
		RelativeInformationScore infScore = new RelativeInformationScore();
		System.out.println(infScore.calculateScore(expected));

		
		//motif = GAATTCG
		//RandomizedGreedyFinder greedy = new RandomizedGreedyFinder(alph, slist, 7, false);
		//Sequence motif = greedy.findMotifs();
		
		GibbsSamplingFinder gibbs = new GibbsSamplingFinder(alph, slist, 7, 1e-3);
		gibbs.getCurrentProfile().printPwm();
		Sequence motif = gibbs.findMotifs();
		gibbs.runMultiple(10);
		
//		RandomProjectionFinder pro = new RandomProjectionFinder(alph, slist, 7, 4, 3, 100);
//		Sequence motif = pro.findMotifs();
//		System.out.println(motif);
	
	}

}
