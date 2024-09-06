/**
 * This class holds all the logic and procedures for hitting a monkey in the face.
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

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Random;

import ilog.concert.*;
import ilog.cplex.*;
import ilog.cplex.IloCplex.UnknownObjectException;

import java.util.ArrayDeque;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Queue;

public class AlgorithmHandler {


	int[] xstar;						// Optimal solution
	int[] probeStar;					// Optimal Probe solution
	
	IloCplex cplex;						// Cplex model
	IloCplex cplexB;					// Cplex model

	// For 2 stage problem
	double LB;
	double UB;

	// For the 2 phase
	int[] probed; //1 if probed in the first stage
	boolean stop;

	int timeLimit = 3600;

	public AlgorithmHandler(DataHandler data) throws InterruptedException, IloException {
		LB = 0;
		UB = 999999;
		xstar = new int[data.n];
		cplex = new IloCplex();
		cplexB = new IloCplex();
		probeStar = new int[data.n];

	}

	public void solveNoProbing(DataHandler data) throws IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, timeLimit);
		cplex.setOut(null);

		IloNumVar[] x = new IloNumVar[data.n];   	
		
		//Create variables 
		for (int i = 0; i < data.n; i++) {
			x[i]= cplex.numVar(0, 1, IloNumVarType.Bool);		
		}

		// KP constraint
		IloLinearNumExpr expr0 = cplex.linearNumExpr();			
		for (int j = 0; j <data.n; j++) {
			expr0.addTerm(data.weight[j], x[j]);
		}
		cplex.addLe(expr0, data.rhs);

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for (int i = 0; i < data.n; i++) {
			obj.addTerm((1-data.probFailure[i])*data.reward[i], x[i]);			
		}

		cplex.addMaximize(obj,"Reward");
		cplex.solve();


		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		No Probing is infeasible!");
		}
		else {
			if(cplex.getStatus() == IloCplex.Status.Feasible || cplex.getStatus() == IloCplex.Status.Optimal){
				data.noProbingObj = data.round(cplex.getObjValue());
				System.out.println("No Probing obj val: "+ data.noProbingObj);
			}
		}

		cplex.clearModel();
		cplex.end();

	}

	public void solvePerfectInfo(DataHandler data) throws InterruptedException, IloException {
		// Bilevel handler
		AlgorithmHandlerBilevel alg = new AlgorithmHandlerBilevel(data);
		
		//Solve perfect info for the VIs
		alg.solvePerfectInfo(data, this);
		
		data.perfectInfoObj = alg.perfectInfoSol.estimatedReward;
	
	}

	public void bilevelValueFunction(DataHandler data, double atime) throws InterruptedException, IloException {

		AlgorithmHandlerBilevel alg = new AlgorithmHandlerBilevel(data);

		//iterarion 0
		int iter=0;
		int nsol = 0;
		// Sample the HPP: compute the solution with no probing to start the sample
		alg.sampleHPPCPLEX(data, this);
		// Update the incumbent solution
		alg.getIncumbentSol(data, probeStar);

		//Start the main loop
		while(alg.UB - alg.LB >alg.tol && (System.currentTimeMillis()-atime)/1000<timeLimit){	
			// Get sample size
			nsol = alg.sample.size();
			//Solve REHPP
			alg.REHPP(data, this, atime);
			// Solve follower's subproblem
			if((System.currentTimeMillis()-atime)/1000<timeLimit) {alg.solveFollower(data,nsol, this);}
			// Update LB
			if(alg.sample.get(nsol).actualReward>alg.LB){ 
				alg.LB = alg.sample.get(nsol).actualReward; 
			};
			System.out.println("Iter: "+iter+" UB: "+alg.UB+" LB: "+alg.LB+" Sample size: "+nsol);
			iter++;
		}	

		// Print the results
		alg.getIncumbentSol(data, probeStar);
		LB = alg.LB;
		UB = alg.UB;
		data.bilevelNumCuts = nsol;

	}

	public void bilevelBB(DataHandler data, double atime) throws InterruptedException, IloException {

		AlgorithmHandlerBilevel alg = new AlgorithmHandlerBilevel(data);

		double[] meanProbe = alg.initializeBB(data, atime, this);		//Solve root node
		if((System.currentTimeMillis()-atime)/1000<timeLimit) { alg.primalHeuristic(data, atime, this, meanProbe);}	//Run primal heuristic
		if((System.currentTimeMillis()-atime)/1000<timeLimit) { alg.runBB(data, atime, this);}	// Run Branch and bound
		LB = alg.LB;
		UB = alg.UB;
		data.bilevelNumCuts = alg.sample.size();
		data.BBNumNodes = alg.numNodesExplored;


	}

	public void runHeuristics(DataHandler data, double atime) throws InterruptedException, IloException {
		// Bilevel handler
		AlgorithmHandlerBilevel alg = new AlgorithmHandlerBilevel(data);
		atime = System.currentTimeMillis();
		
//		alg.greedyHeuristic(data, atime, this);
		data.heuTime1 = ((System.currentTimeMillis()-atime)/1000.0);
		
		alg.greedyHeuristicPI(data, atime, this);
		data.heuTime2 = ((System.currentTimeMillis()-atime)/1000.0)-data.heuTime1;
		
//		alg.heuristicER2(data, atime, this);
		data.heuTime3 = ((System.currentTimeMillis()-atime)/1000.0)-(data.heuTime1+data.heuTime2);	
	}

	// Auxiliary class for sorting
	class Item implements Comparable<Item>
	{
		int pos;
		double val;
	
		public Item(int pos, double val) {
			this.pos = pos;
			this.val = val;
		}
		@Override
		public int compareTo(Item o) {
			return Double.compare(val, o.val);
		}
	
	}

	public void ComputeBounds(DataHandler data) {

		for (int i = 0; i < data.n; i++) {
			data.boundPV += 0.25*data.reward[i];
		}
		
		System.out.println("BOUND ON PROBING VALUE: "+data.boundPV);
		
		
		for (int i = 0; i < data.n; i++) {
			if(probeStar[i]==0) {data.boundPG += 0.5*data.reward[i];}
			else {System.out.println("WEPA "+i);}
		}
		
		System.out.println("BOUND ON PRICE GAP: "+data.boundPG);
		
		
	}


}
