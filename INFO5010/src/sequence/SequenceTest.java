package sequence;

import org.junit.Test;

public class SequenceTest
{

	@Test
	public void test()
	{
		Alphabet alph = new Alphabet("ACGT","");
		Sequence seq = Sequence.generateRandomSequence(alph, 30);
		System.out.println("-Original-");
		System.out.println(seq);
		System.out.println("-Mutation 1-");
		seq.mutate(10);
		System.out.println(seq);
		System.out.println("-Mutation 2-");
		seq.mutate(0.1);
		System.out.println(seq);
		
	}

}
