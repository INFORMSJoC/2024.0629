/**
 * This is the main class.
 * 
 * Ref.: TBD
 * 
 * @author L. Lozano & D. Bergman
 * @url http://business.uc.edu/academics/departments/obais/faculty/leonardo-lozano.html
 * 
 */
package KP;

import ilog.concert.IloException;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;

public class Main {

	public static void main(String[] args) throws IOException, InterruptedException, IloException {
		java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("Results.txt", true));
		// BILEVEL
		int[] S = {10, 30, 50}; //Scenarios
		int[] budget = {1, 3, 6}; //Probing budgets
		for (int q = 0; q < S.length; q++) {
			for (int b = 0; b < budget.length; b++) {
				for (int s = 0; s < 10; s++) {
					// Begin the time count						
					double Atime = System.currentTimeMillis();

					DataHandler data = new DataHandler();
					data.genInstance(20, s, 0.1); // Create random instance: n, seed, rhs mult
					//data.outputInstance(s);
					data.genScenarios(S[q], budget[b]); // Generate scenarios: number of scenarios, probing budget
					//data.genAllScenarios(budget[b]); // Generate scenarios: number of scenarios, probing budget

					// Create an AlgorithmHandler
					AlgorithmHandler alg = new AlgorithmHandler(data);
					alg.solveNoProbing(data);
					alg.solvePerfectInfo(data);
					alg.bilevelValueFunction(data, Atime);
					//alg.runHeuristics(data, Atime);
					

					System.out.println("       EXECUTION TIME: "+((System.currentTimeMillis()-Atime)/1000.0)+" "+alg.LB+" "+alg.UB);
					ps.println("VF "+(data.n)+" "+data.numScenarios+" "+data.budget+" "+s+" "+((System.currentTimeMillis()-Atime)/1000.0)+" "+alg.LB+" "+alg.UB+" "+data.bilevelNumCuts+" "+data.noProbingObj+" "+data.perfectInfoObj);
					//ps.println("Heuristics "+(data.n)+" "+data.numScenarios+" "+data.budget+" "+s+" "+data.heuTime1+" "+data.heuTime2+" "+data.heuTime3+" "+data.noProbingObj+" "+data.perfectInfoObj+" "+data.greedyObj1+" "+data.greedyObjPI+" "+data.heuristicER+" "+data.budgetUB);
					alg = null;
					System.gc();
				}
			}
		}
	}



	}
