package phm;

import java.io.IOException;

/**
 * Example of how to use the PHM algorithm from the source code.
 * @author Philippe Fournier-Viger, 2016
 */
public class Run {

	public static void main(String [] arg) throws IOException{

		String output = "output";

		// =======================
		// EXAMPLE FROM THE ARTICLE : 
		String input = "input";
		int min_utility = 20;   
		int minPeriodicity = 1;  // minimum periodicity parameter (a number of transactions)
		int maxPeriodicity = 3;  // maximum periodicity parameter (a number of transactions)
		int minAveragePeriodicity = 1;  // minimum average periodicity (a number of transactions)
		int maxAveragePeriodicity = 2;  // maximum average periodicity (a number of transactions)
		// =======================

		// Applying the PHM algorithm
		AlgoPHM algorithm = new AlgoPHM();
		// To disable some optimizations:
		//algorithm.setEnableEUCP(false); 
		//algorithm.setEnableESCP(false);
		
		// Run the algorithm
		algorithm.runAlgorithm(input, output, min_utility, 
				minPeriodicity, maxPeriodicity, minAveragePeriodicity, 
				maxAveragePeriodicity);
		
		algorithm.printStats();
		
	}

}
