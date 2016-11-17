package fhm;

import java.io.IOException;

/**
 * Example of how to use the FHM algorithm 
 * from the source code.
 * @author Philippe Fournier-Viger, 2014
 */
public class Run {

	public static void main(String [] arg) throws IOException{
		
		String input = "input";
		String output = "output";

		int min_utility = 30;  // 
		
		// Applying the HUIMiner algorithm
		AlgoFHM fhm = new AlgoFHM();
		fhm.runAlgorithm(input, output, min_utility);
		fhm.printStats();

	}
}
