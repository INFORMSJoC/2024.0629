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
package KP;

import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.concert.IloNumVarType;
import ilog.cplex.IloCplex;


public class DataHandler {

	// Number of nodes
	public int n;

	//Rewards
	int [] reward;
	// Weights 
	int [] weight;
	int rhs;

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
	double heuTime1;
	double heuTime2;
	double heuTime3;
	int budgetUB;
	
	double boundPV;
	double boundPG;

	// Create empty datahandler
	public DataHandler() {
		BBNumNodes = 0;
	}

	public void genInstance(int _n, int seed, double alpha) throws IloException {
		r = new Random(seed);
		n = _n; // include depot
		weight = new int[n];
		reward = new int[n];

		int sum = 0;
		for (int i = 0; i < n; i++) {
			weight[i] = 1+r.nextInt(50);
			sum += weight[i];
			reward[i] = 50+r.nextInt(51);
		}
		rhs = (int) (sum*alpha);
		System.out.println("INSTANCE: "+Arrays.toString(weight)+" <= "+rhs);

	}

	public void genScenarios(int numS, int _budget) {
		// Generate scenarios
		numScenarios = numS;
		budget = _budget;
		scenarios = new int[numScenarios][n];
		probFailure = new double[n];
		for (int k = 0; k < numScenarios; k++) {
			for (int i = 0; i < n; i++) {
				if(r.nextDouble() <= 0.5) {
					scenarios[k][i] = 1;
					probFailure[i]++;
				}
			}
			//System.out.println("Scenario "+k+": "+Arrays.toString(scenarios[k]));
		}

		//Compute probabilities
		for (int i = 0; i < n; i++) {
			probFailure[i] = probFailure[i]/numScenarios;
			//System.out.println("Probability jammer at "+i+": "+probFailure[i]+" reward: "+reward[i]+" estimated: "+(1-probFailure[i])*reward[i]);
		}
	}

	public void genAllScenarios(int _budget) {
		// Generate scenarios
		numScenarios = (int) Math.pow(2, n);
		budget = _budget;
		scenarios = new int[numScenarios][n];
		probFailure = new double[n];
		ArrayList<Integer> aux = new ArrayList<Integer>();
		aux.add(0);
		recursionAllScenarios(new int[n], 0, aux);

		//Compute probabilities
		for (int k = 0; k < numScenarios; k++) {
			for (int i = 0; i < n; i++) {
				if(scenarios[k][i] == 1) {
					probFailure[i]++;
				}
			}
			//System.out.println("Scenario "+k+": "+Arrays.toString(scenarios[k]));
		}
		for (int i = 0; i < n; i++) {
			probFailure[i] = probFailure[i]/numScenarios;
			System.out.println("Probability jammer at "+i+": "+probFailure[i]+" reward: "+reward[i]+" estimated: "+(1-probFailure[i])*reward[i]);
		}
	}
	
	public void recursionAllScenarios(int[] s, int i, ArrayList<Integer> count)
	{
		if (i == n) 
		{
			for (int j = 0; j < n; j++) {
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

	public void outputInstance(int s) throws FileNotFoundException {
		java.io.PrintStream ps = new java.io.PrintStream( new java.io.FileOutputStream("KP_"+n+"_"+s+".txt"));
		ps.println("profit_coefficients:");
		for (int i = 0; i < n; i++) {
			ps.print(reward[i]+" ");
		}
		ps.println("");
		ps.println("investment_coefficients:");
		for (int i = 0; i < n; i++) {
			ps.print(weight[i]+" ");
		}
		ps.println("");
		ps.println("total_funds:");
		ps.println(rhs);

		
	}




}

