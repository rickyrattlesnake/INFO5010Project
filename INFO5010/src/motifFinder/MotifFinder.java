package motifFinder;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Scanner;

import motifsearch.Finder;
import motifsearch.GibbsSamplingFinder;
import motifsearch.RandomProjectionFinder;
import motifsearch.RandomizedGreedyFinder;
import scoring.ExpectationScore;
import scoring.ExpectedInformationScore;
import scoring.FrequencyScore;
import scoring.RelativeInformationScore;
import scoring.Score;
import sequence.Alphabet;
import sequence.Profile;
import sequence.Sequence;

public class MotifFinder
{
	private int motifLength = 0;
	private Alphabet alphabet = new Alphabet("ACGT",""); 				//Default DNA alphabet
	private List<Sequence> seqList = new ArrayList<Sequence>();			//Default empty list of sequences
	private int numTrials = 1;											//Default 1 trial
	private Sequence consensusMotif = null;
	private Map<Sequence, Integer> alignments = null; 					//vector of motif start positions for each sequence
	private Profile profile = null;
	private Finder algorithm = null;									//chosen motif Finding algorithm
	private Score scorer = null; 										//chosen scoring metric, if null use default
	private Map<Sequence, Integer> perfectAlignments = null; 			//memory of inserted motif alignments
	
	private static String HELP_TEXT = 
			  "*-Help Text-*\n"
			+ "\n"
			+ "-- Commands --\n"
			+ "\n"
			+ "load-file <inputFile> 							: loads the alphabet, symbol distribution, motif length and sequences from a file\n"
			+ "generate <length> <quantity>						: generates the specified quantity of random sequences of specified length (assumes alphabet = 'ACGT')\n"
			+ "insert-motif <length> <mutation-rate> 			: inserts motifs into sequences with a given mutation-rate (%) per base and length\n"
			+ "trials <number> 									: number of trials to run for each finder [default = 1]\n"
			+ "clear 											: clears all motifFinder parameters and sequence lists\n"
			+ "\n"
			+ "print-motif 								: prints the consensus motif\n"
			+ "print-profile 							: prints the current profile probability matrix\n"
			+ "print-alignments 						: prints the current alignment vector\n"
			+ "print-sequences 							: prints all sequences\n"
			+ "print-inserted-motif						: prints insert-consensus motif if inserted\n"
			+ "\n"
			+ "find-motif <algorithm> <param1> <param2> ... : runs the given algorithm for given params. Check Below for details\n"
			+ "set-scoring <scoringType> [param1]   		: sets the scoring metric with an optional param. Check Below for details\n"
			+ "\n"
			+ "\n-- Algorithms and Required Parameters --\n"
			+ "greedy <updateStep>											: Runs the randomized greedy finder, updateStep=true|false if updates are to happen at each alignment change. Recommended this be set to true\n"
			+ "gibbs <optimizationThreshold>								: Runs the gibbs sampling finder, optimisationThresh determines when to stop optimisation (recommend = 1e-7) \n"
			+ "projection <projectionSize> <binThreshold> <numIterations>	: Runs the random projection finder; projectionSize refers to the size of the hashed kmer, binThreshold determines which bins are selected for further anlysis, numIterations determines how many k-l templates are projected\n"
			+ "\n-- Types of Scoring Metrics--\n"
			+ "frequency 							: simple frequency summation to measure the strength of consensus\n"
			+ "expectation [pseudoZero] 			: sum(log_2(p)) using the position probability matrix\n"
			+ "expected-information [pseudoZero] 	: sum(p * log_2(p)) using the position probability matrix\n"
			+ "relative-information [pseudoZero] 	: sum(p * log_2(p/b)) using the position weight matrix\n"
			+ "";
	//TODO: fill in algorithms and parameters
			
	public static void main(String[] args)
	{
		MotifFinder mFinder = new MotifFinder();
		Scanner sc = new Scanner(System.in);
		
		//Print out a help indicator on first run
		System.out.println("Type 'help' for help text");
		
		while(true)
		{
			System.out.print("Enter Command >> ");
			String[] input = sc.nextLine().trim().split(" ");
			String command = input[0];
			switch(command)
			{
			case "help":
				System.out.println(HELP_TEXT);
				break;
			case "load-file":
				if(input.length == 2)
				{
					mFinder.loadFile(input[1]);
				}else
				{
					System.err.println("Incorrect arguments : type 'help' for help-text");
				}
				break;
			case "generate":
				if(input.length == 3)
				{
					try{
						mFinder.generate(Integer.parseInt(input[1]), Integer.parseInt(input[2]));
					}catch(NumberFormatException ex)
					{
						System.err.println("Incorrect arguments : type 'help' for help-text");
					}
				}else
				{
					System.err.println("Incorrect arguments : type 'help' for help-text");
				}
				break;
			case "insert-motif":
				if(input.length == 3)
				{
					try{
						mFinder.insertMotif(Integer.parseInt(input[1]), Double.parseDouble(input[2]));
					}catch(NumberFormatException | ArrayIndexOutOfBoundsException ex)
					{
						System.err.println("Incorrect arguments : type 'help' for help-text");
					}
				}else
				{
					System.err.println("Incorrect arguments : type 'help' for help-text");
				}
				break;
			case "trials":
				if(input.length == 2)
				{
					try{
						mFinder.setNumTrials(Integer.parseInt(input[1]));
					}catch(NumberFormatException ex)
					{
						System.err.println("Incorrect arguments : type 'help' for help-text");
					}
				}else
				{
					System.err.println("Incorrect arguments : type 'help' for help-text");
				}
				break;
			case "clear":
				mFinder.clear();
				break;
			case "print-motif":
				mFinder.printMotif();
				break;
			case "print-profile":
				mFinder.printProfile();
				break;
			case "print-alignments":
				mFinder.printAlignments();
				break;
			case "print-sequences":
				mFinder.printAllSequences();
				break;
			case "print-inserted-motif":
				mFinder.printInsertedMotif();
				break;
			case "find-motif":
				try{
					mFinder.runFinder(input);
				}catch(NumberFormatException ex)
				{
					System.err.println("Incorrect arguments : type 'help' for help-text");
				}
				break;
			case "set-scoring":
				try{
					mFinder.setScorer(input);
				}catch(NumberFormatException ex)
				{
					System.err.println("Incorrect arguments : type 'help' for help-text");
				}
				break;
			default:
				System.err.println("Command not recognised : type 'help' for help text");				
			}
		}
	}
	
	
	/**
	 * Generates the specified quantity of random sequences of specified length. 
	 * Assumes alphabet = 'ACGT' and a uniform background probability for each DNA base
	 * @param quantity
	 * @param length
	 */
	public void generate(int length, int quantity)
	{
		Alphabet alph = new Alphabet("ACGT", "");
		for(int q=0; q < quantity; ++q)
		{
			Sequence randSeq = Sequence.generateRandomSequence(alph, length);
			seqList.add(randSeq);
		}
	}
	
	
	/**
	 * Inserts motifs into all existing sequences with a given mutation-rate per base.
	 * Returns a list of starting positions for each motif and corresponding sequence.
	 * Assigned the perfect alignment map.
	 * @param mutationRate
	 * @return startingPos
	 */
	public void insertMotif(int length, double mutationRate)
	{
		motifLength = length; 
		Alphabet alph = new Alphabet("ACGT", "");
		int[] startPos = new int[seqList.size()];
		
		
		Sequence motif = Sequence.generateRandomSequence(alph, length);
		for(int i=0; i < seqList.size(); ++i)
		{
			Sequence mutatedMotif = motif.copy();
			mutatedMotif.mutate(mutationRate);
			startPos[i] = seqList.get(i).insertMotif(mutatedMotif);
		}
		
		perfectAlignments = new HashMap<Sequence, Integer>();
		System.out.print("Inserted Motif positions: [");
		for(int i=0; i < startPos.length; ++i)
		{
			System.out.print(startPos[i] + ", ");
			perfectAlignments.put(seqList.get(i), startPos[i]);
		}
		System.out.print("]\n");
	}
	
	/**
	 * Clears the sequence list and sets all parameters to default values.
	 * NB: clears the sequence list in place, so any dependent finders will also lose the sequences
	 */
	public void clear()
	{
		alphabet = new Alphabet("ACGT", "");
		motifLength = 0;
		seqList.clear();
		numTrials = 1;
		consensusMotif = null;
		alignments = null;
		profile = null;
		scorer = null;
		perfectAlignments = null;
	}
	
	/**
	 * Creates an instance of the specified finder and runs the specified number of trials
	 * @param input 
	 */ 
	public void runFinder(String[] input)
	{
		if(input.length < 2)
		{
			System.err.println("Incorrect arguments : type 'help' for help-text");
			return;
		}
		if(motifLength == 0 || seqList.isEmpty())
		{
			System.err.println("Sequences and Motif lengths must be initialized first");
			return;
		}
		
		switch(input[1])	//second argument is the algorithm-type
		{
		case "greedy":
			if(input.length != 3)
			{
				System.err.println("Incorrect arguments : type 'help' for help-text");
				return;
			}
			if(scorer == null)
				algorithm = new RandomizedGreedyFinder(alphabet, seqList, motifLength, Boolean.parseBoolean(input[2]));
			else
				algorithm = new RandomizedGreedyFinder(alphabet, seqList, motifLength, Boolean.parseBoolean(input[2]), scorer);
			break;
		case "gibbs":
			if(input.length != 3)
			{
				System.err.println("Incorrect arguments : type 'help' for help-text");
				return;
			}
			if(scorer == null)
				algorithm = new GibbsSamplingFinder(alphabet, seqList, motifLength, Double.parseDouble(input[2]));
			else
				algorithm = new GibbsSamplingFinder(alphabet, seqList, motifLength, Double.parseDouble(input[2]), scorer);
			break;
		case "projection":
			if(input.length != 5)
			{
				System.err.println("Incorrect arguments : type 'help' for help-text");
				return;
			}
			if(scorer == null)
				algorithm = new RandomProjectionFinder(alphabet, seqList, motifLength, Integer.parseInt(input[2]), Integer.parseInt(input[3]), Integer.parseInt(input[4]));
			else
				algorithm = new RandomProjectionFinder(alphabet, seqList, motifLength, Integer.parseInt(input[2]), Integer.parseInt(input[3]), Integer.parseInt(input[4]), scorer);
			break;
		default:
			System.err.println("Incorrect Algorithm Type : type 'help' for help-text");
			return;
		}
		
		
		//Run the number of trials for the finder
		consensusMotif = algorithm.runMultiple(numTrials);
		//Set the profile, and alignments
		profile = algorithm.getCurrentProfile();
		alignments = profile.getAlignmentStarts();
	}
	
	/**
	 * prints the consensus motif if one has been set
	 */
	public void printMotif()
	{
		if(consensusMotif == null)
		{
			System.err.println("Motif was not set, ensure you have loaded the data and parameters and run a finder.");
		}
		else
		{
			System.out.println(consensusMotif);
		}
	}
	
	/**
	 * prints the consensus motif if one has been set
	 */
	public void printInsertedMotif()
	{
		if(perfectAlignments == null)
		{
			System.err.println("Motif was not inserted.");
		}
		else
		{
			Profile prof = new Profile(alphabet, seqList, motifLength, perfectAlignments);
			System.out.println(prof.getConsensus());
		}
	}
	
	/**
	 * Prints all sequences. One per line.
	 */
	public void printAllSequences()
	{
		if(seqList == null)
		{
			System.err.println("No sequences have been loaded.");
		}
		else
		{
			for(Sequence seq : seqList)
				System.out.println(seq);
		}
	}
	
	/**
	 * Prints the alignment vector, assumes that the sequence list has
	 * not been modified since the last run of Finder
	 */
	public void printAlignments()
	{
		if(alignments == null)
		{
			System.err.println("Alignments were not set, ensure you have loaded the data and parameters and run a finder.");
		}
		else
		{
			for(Sequence seq : seqList)
				System.out.println(alignments.get(seq));
		}
	}
	
	/**
	 * Prints the alignment vector, assumes that the sequence list has
	 * not been modified since the last run of Finder
	 */
	public void printProfile()
	{
		if(profile == null)
		{
			System.err.println("Profile was not set, ensure you have loaded the data and parameters and run a finder.");
		}
		else
		{
			System.out.println(profile);
		}
	}
	
	/**
	 * Loads alphabet, motif length and sequences from a file
	 * 
	 * File format: 
	 * 
	 * <Alphabet> 									e.g. ACGT
	 * <probability distribution for each symbol>	e.g. 0.25 0.3 0.25 0.20
	 * <Motif Length>								e.g. 6
	 * <sequence1>									e.g. ACGATATCTTCGTAGCTAG
	 * <sequence2>									e.g. ACTGATCGATGCTAGATCG
	 * ...
	 * 
	 * NB: no newlines before or at the end of the expected input 
	 * 
	 * @param fileName	input-file name, assumes current directory
	 * @param seqList	list to load sequences from the file in order of appearance
	 * @return motif-length length of the motif specified 
	 */
	public void loadFile(String fileName)
	{
		//TODO: error handling
		BufferedReader br = null;
		Scanner sc = null;
		try
		{
			br = new BufferedReader(new FileReader(fileName));
			sc = new Scanner(br);
		} 
		catch (FileNotFoundException e)
		{
			e.printStackTrace();
			System.exit(1);
		}
		//Extract the Alphabet
		String alphString = sc.nextLine().trim();
		//Extract the symbol prob. distribution
		double[] probDist = new double[alphString.length()];
		String[] probStrings = sc.nextLine().trim().split(" ");
		for(int i=0; i < probDist.length; ++i)
		{
			probDist[i] = (Double.parseDouble(probStrings[i]));
		}
		//Extract the motif-length
		motifLength = Integer.parseInt(sc.nextLine().trim());
		//Create an alphabet 
		alphabet = new Alphabet(alphString, "", probDist);
		//Create a Sequence list from each sequence
		while(sc.hasNext())
		{
			String line = sc.nextLine().trim();
			if(line == "")
				continue;
			
			String seqData = line;
			Sequence seq = new Sequence(alphabet, seqData);
			seqList.add(seq);
		}
	}
	
	public void setNumTrials(int trials)
	{
		numTrials = trials;
	}
	
	/**
	 * Determines which scorer to set with optional parameters
	 * @param input
	 */
	public void setScorer(String[] input)
	{
		if(input.length == 2)
		{
			switch(input[1])	//second argument is the score-type
			{
			case "frequency":
				scorer = new FrequencyScore();
				break;
			case "relative-information":
				scorer = new RelativeInformationScore();
				break;
			case "expected-information":
				scorer = new ExpectedInformationScore();
				break;
			case "expectation":
				scorer = new ExpectationScore();
				break;
			default:
				System.err.println("Incorrect Scoring Type : type 'help' for help-text");
			}
		}
		else if(input.length == 3)
		{
			double pseudoZero = Double.parseDouble(input[2]);
			switch(input[1])	//second argument is the score-type
			{
			case "relative-information":
				scorer = new RelativeInformationScore(pseudoZero);
				break;
			case "expected-information":
				scorer = new ExpectedInformationScore(pseudoZero);
				break;
			case "expectation":
				scorer = new ExpectationScore(pseudoZero);
				break;
			default:
				System.err.println("Incorrect Scoring Type : type 'help' for help-text");
			}
		}
		else
		{
			System.err.println("Incorrect arguments : type 'help' for help-text");
		}
	}
}
