/**
 * This is the main class.
 * 
 * Ref.: TBD
 * 
 * @author L. Lozano & D. Bergman
 * @url http://business.uc.edu/academics/departments/obais/faculty/leonardo-lozano.html
 * 
 */
package SP;

import ilog.concert.IloException;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

public class Main {

	public static void main(String[] args) throws IOException, InterruptedException, IloException {
		java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("Results.txt", true));
		// BILEVEL
		int[] S = {10, 30, 50};
		int[] budget = {4, 8, 16};
		for (int q = 0; q < S.length; q++) {
			for (int b = 0; b < budget.length; b++) {
				for	(int s = 0; s < 10 ; s++) {
					// Begin the time count						
					double Atime = System.currentTimeMillis();

					DataHandler data = new DataHandler();
					data.ReadDimacs("data/Grid5x5_100-100_"+s+".txt");
					data.genScenarios(S[q], budget[b], s); // Generate scenarios: number of scenarios, probing budget
					//data.genAllScenarios(budget[b]); // Generate scenarios: number of scenarios, probing budget

					// Create an AlgorithmHandler
					AlgorithmHandler alg = new AlgorithmHandler(data);
					alg.solveNoProbing(data);
					alg.solvePerfectInfo(data);
					alg.bilevelValueFunction(data, Atime);
					//alg.runHeuristics(data, Atime);
					

					System.out.println("       EXECUTION TIME: "+((System.currentTimeMillis()-Atime)/1000.0)+" "+alg.LB+" "+alg.UB);
					ps.println("VF "+(data.CsvInput)+" "+data.numScenarios+" "+data.budget+" "+((System.currentTimeMillis()-Atime)/1000.0)+" "+alg.LB+" "+alg.UB+" "+data.bilevelNumCuts+" "+data.noProbingObj+" "+data.perfectInfoObj);
					//ps.println("Heuristics "+(data.CsvInput)+" "+data.numScenarios+" "+data.budget+" "+data.heuTime1+" "+data.heuTime2+" "+data.heuTime3+" "+data.noProbingObj+" "+data.perfectInfoObj+" "+data.greedyObj1+" "+data.greedyObjPI+" "+data.heuristicER+" "+data.budgetUB);
					alg = null;
					System.gc();

				}
			}
		}
	}



}
