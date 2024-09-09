/**
 * This class holds all the data.
 * 
 * Ref.: Lozano, L. and Smith, J. C. (2015). 
 * A Sampling-Based Exact Approach for the Bilevel Mixed Integer Programming Problem
 * 
 * @author L. Lozano & J. C. Smith
 * @affiliation Clemson University
 * @url www.leo-loza.com
 * 
 */
package SP;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.StringTokenizer;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.cplex.IloCplex;

public class DataHandler {

	// Name of the instance
	String CsvInput;
	// Number of arcs
	int numArcs;
	// Number of nodes
	public static int numNodes;
	// Destinantion node
	public int lastNode;
	// Source node
	public int source;
	// All the arcs in the network stored in a vector where Arcs[i][0]= Tail for arc i and Arcs[i][1]= Head for arc i 
	public static int[][] arcs;
	// The distance attribute for any arc i
	public static int[] cost;
	// The delay attribute for any arc i
	public static int[] delay;
	
	//Scenarios for 2-stage
	int numScenarios;
	int [][] scenarios;
	int budget;
	Random r;
	double[] probFailure;
	
	// Performance measures
	int BBNumNodes;
	int bilevelNumCuts;
	double noProbingObj;
	double perfectInfoObj;
	double greedyObj1;
	double greedyObjPI;
	double heuristicER;
	double heuristicEV;
	double heuTime1;
	double heuTime2;
	double heuTime3;
	int budgetUB;
	double boundVI;
	double boundPG;
	
	// Create empty datahandler
	public DataHandler() {

	}

	// This procedure reads data from a data file in DIMACS format
	public void ReadDimacs(String input) throws NumberFormatException, IOException {
		CsvInput = input;
		File file = new File(CsvInput);

		BufferedReader bufRdr = new BufferedReader(new FileReader(file));
		// Read first line
		String line = bufRdr.readLine();
		String[] array = line.split(" ");
		numNodes = Integer.parseInt(array[2]);
		numArcs = Integer.parseInt(array[3]);
		System.out.println("Network size: "+numNodes+" "+numArcs);
		
		// Read Second line
		line = bufRdr.readLine();
		array = line.split(" ");
		source = Integer.parseInt(array[1]);

		// Read third line
		line = bufRdr.readLine();
		array = line.split(" ");
		lastNode = Integer.parseInt(array[1]);

		//System.out.println("WEPAAAA: "+source+" "+lastNode);
		
		arcs = new int[numArcs][2];
		cost = new int[numArcs];
		delay = new int[numArcs];
				
		String[] readed = new String[5];

		int row = 0;
		int col = 0;
		
		while ((line = bufRdr.readLine()) != null && row < numArcs + 3) {
			StringTokenizer st = new StringTokenizer(line, " ");
			while (st.hasMoreTokens()) {
				// get next token and store it in the array
				readed[col] = st.nextToken();
				col++;
			}

			if (row >= 0) {
				arcs[row][0] = (Integer.parseInt(readed[1]) );
				arcs[row][1] = (Integer.parseInt(readed[2]) );
				cost[row] = Integer.parseInt(readed[3]);
				delay[row] = Integer.parseInt(readed[4]);
			}

			col = 0;
			row++;

		}
		
	}


	public void genScenarios(int numS, int _budget, int seed) {
		// Generate scenarios
		r = new Random(seed);
		numScenarios = numS;
		budget = _budget;
		scenarios = new int[numScenarios][numArcs];
		probFailure = new double[numArcs];
		for (int k = 0; k < numScenarios; k++) {
			for (int i = 0; i < numArcs; i++) {
				if(r.nextDouble() <= 0.5) {
					scenarios[k][i] = 1;
					probFailure[i]++;
				}
			}
			//System.out.println("Scenario "+k+": "+Arrays.toString(scenarios[k]));
		}
		
		//Compute probabilities
		for (int i = 0; i < numArcs; i++) {
			probFailure[i] = probFailure[i]/numScenarios;
			//System.out.println("Probability jammer at "+i+": "+probFailure[i]+" reward: "+cost[i]+" estimated: "+(1-probFailure[i])*cost[i]);
		}
		
		
	}
	
	public void genAllScenarios(int _budget) {
		// Generate scenarios
		numScenarios = (int) Math.pow(2, numArcs);
		budget = _budget;
		scenarios = new int[numScenarios][numArcs];
		probFailure = new double[numArcs];
		ArrayList<Integer> aux = new ArrayList<Integer>();
		aux.add(0);
		recursionAllScenarios(new int[numArcs], 0, aux);

		//Compute probabilities
		for (int k = 0; k < numScenarios; k++) {
			for (int i = 0; i < numArcs; i++) {
				if(scenarios[k][i] == 1) {
					probFailure[i]++;
				}
			}
			//System.out.println("Scenario "+k+": "+Arrays.toString(scenarios[k]));
		}
		for (int i = 0; i < numArcs; i++) {
			probFailure[i] = probFailure[i]/numScenarios;
			System.out.println("Probability jammer at "+i+": "+probFailure[i]+" cost: "+cost[i]+" estimated: "+(1-probFailure[i])*cost[i]);
		}
	}

	public void recursionAllScenarios(int[] s, int i, ArrayList<Integer> count)
	{
		if (i == numArcs) 
		{
			for (int j = 0; j < numArcs; j++) {
				scenarios[count.get(0)][j] = s[j];
			}
			count.set(0, count.get(0)+1);
			//System.out.println(Arrays.toString(s));
			return;
		}

		// First assign "0" at ith position
		// and try for all other permutations
		// for remaining positions
		s[i] = 0;
		recursionAllScenarios(s, i + 1, count);

		// And then assign "1" at ith position
		// and try for all other permutations
		// for remaining positions
		s[i] = 1;
		recursionAllScenarios(s, i + 1, count);
	}
	public double round(double value) {
		double rounded;
		rounded = Math.round(value*10000)/10000.0;
		return rounded;
	}




}

