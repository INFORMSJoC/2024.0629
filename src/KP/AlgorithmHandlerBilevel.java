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


import ilog.concert.*;
import ilog.cplex.*;
import ilog.cplex.IloCplex.UnknownObjectException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Random;


public class AlgorithmHandlerBilevel {


	double LB; 								// Best lower bound
	double UB; 								// Best upper bound

	double tol;								// Tolerance
	double M;								// Big M
	ArrayList<SolutionBilevel> sample;		// Current sample
	SolutionBilevel noProbeSol;				// Solution when there is no probing
	SolutionBilevel perfectInfoSol;			// Solution when there is perfect info
	SolutionBilevel xstar;					// Optimal sol

	IloCplex cplex;							// Cplex model 1
	IloCplex cplex2;						// Cplex model 2

	int maxVisits;

	ArrayList<BBNode> tree;
	int numNodesExplored;

	public AlgorithmHandlerBilevel(DataHandler data) throws InterruptedException, IloException {
		UB = 999999999;
		LB = 0;
		sample = new ArrayList<SolutionBilevel>();
		tol = 0.001;
		xstar = new SolutionBilevel(data);
		maxVisits = 999999;
		cplex = new IloCplex();
		cplex2 = new IloCplex();
		numNodesExplored = 0;
	}



	public int vectorMult(int[] v, int[] x) {
		int result = 0;

		for (int i = 0; i < x.length; i++) {
			result+=v[i]*x[i];
		}

		return result;
	}

	public double round(double x){
		double rounded = Math.floor(x*1000000000)/1000000000.0;

		return rounded;
	}



	public void solvePerfectInfo(DataHandler data, AlgorithmHandler a) throws UnknownObjectException, IloException {
		// Initialize solution with all probing and compute big M
		perfectInfoSol = new SolutionBilevel(data);
		perfectInfoSol.actualReward = -1; //Signal the subproblem that we are solving for perfect info 
		for (int i = 0; i < data.n; i++) {
			perfectInfoSol.probed[i] = 1;
		}

		//Solve for each scenario
		for (int k = 0; k < data.numScenarios; k++) {
			solveSubScenario(data, data.scenarios[k], k, perfectInfoSol, a);
		}
		// Estimated and actual are the same with perfect info!!!!
		perfectInfoSol.estimatedReward = perfectInfoSol.estimatedReward/(data.numScenarios+0.0);
		data.perfectInfoObj = perfectInfoSol.estimatedReward;
		System.out.println("Perfect info obj val: "+data.perfectInfoObj);
		//perfectInfoSol.print();

	}



	public void REHPP(DataHandler data, AlgorithmHandler a, double atime) throws IloException {
		cplex = new IloCplex();
		double timeRemaining = a.timeLimit - (System.currentTimeMillis()-atime)/1000;
		cplex.setParam(IloCplex.Param.TimeLimit, timeRemaining);
		//cplex.setParam(IloCplex.IntParam.MIPDisplay, 0);
		//cplex.setOut(null);

		IloNumVar[][] x = new IloNumVar[data.numScenarios][data.n];   	
		IloNumVar[] probeVar = new IloNumVar[data.n];   										// 1 if node is probed
		IloNumVar[][] binaryMultVar = new IloNumVar[data.numScenarios][data.n];   	// x_ki*probe_ki

		//Create variables 
		for (int i = 0; i < data.n; i++) {
			for (int k = 0; k < data.numScenarios; k++) {
				x[k][i]= cplex.numVar(0, 1, IloNumVarType.Bool);		
			}
		}
		for (int i = 0; i < data.n; i++) {
			probeVar[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		for (int i = 0; i < data.n; i++) {
			for (int k = 0; k < data.numScenarios; k++) {
				binaryMultVar[k][i] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
		}

		// KP constraints
		for (int k = 0; k < data.numScenarios; k++) {
			IloLinearNumExpr expr0 = cplex.linearNumExpr();			
			for (int j = 0; j <data.n; j++) {
				expr0.addTerm(data.weight[j], x[k][j]);
			}
			cplex.addLe(expr0, data.rhs);
		}

		// Probing budget
		IloLinearNumExpr exprB = cplex.linearNumExpr();
		for (int i = 0; i < data.n; i++) {
			exprB.addTerm(1, probeVar[i]);		
		}
		cplex.addLe(exprB, data.budget);


		//////////////////////////////////////////////////////////////////////////////////////////////


		// Add obj value constraints: Leader Sol Obj >= Sample Sol Obj 
		for (int j = 0; j < sample.size(); j++) {
			SolutionBilevel sol = sample.get(j);
			for (int k = 0; k < data.numScenarios; k++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				// Leader Sol obj
				double rhs = 0;
				for (int i = 0; i <data.n; i++) {
					//RHS
					expr.addTerm(data.reward[i]*(1-data.scenarios[k][i]), binaryMultVar[k][i]);
					expr.addTerm(data.reward[i]*(1-data.probFailure[i]), x[k][i]);
					expr.addTerm(-data.reward[i]*(1-data.probFailure[i]), binaryMultVar[k][i]);
					//LHS
					expr.addTerm(-data.reward[i]*(1-data.scenarios[k][i])*sol.xhatf[k][i], probeVar[i]);
					expr.addTerm(data.reward[i]*(1-data.probFailure[i])*sol.xhatf[k][i], probeVar[i]);
					rhs += data.reward[i]*(1-data.probFailure[i])*sol.xhatf[k][i];
				}

				cplex.addGe(expr, rhs-tol);
			}
		}

		//Binary multiplication envelope
		for (int i = 0; i < data.n; i++) {
			for (int k = 0; k < data.numScenarios; k++) {
				IloLinearNumExpr exprA = cplex.linearNumExpr();
				exprA.addTerm(1, binaryMultVar[k][i]);
				exprA.addTerm(-1, x[k][i]);
				cplex.addLe(exprA, 0);

				IloLinearNumExpr exprC = cplex.linearNumExpr();
				exprC.addTerm(1, binaryMultVar[k][i]);
				exprC.addTerm(-1, probeVar[i]);
				cplex.addLe(exprC, 0);

				IloLinearNumExpr exprL = cplex.linearNumExpr();
				exprL.addTerm(1, binaryMultVar[k][i]);
				exprL.addTerm(-1, probeVar[i]);
				exprL.addTerm(-1, x[k][i]);
				cplex.addGe(exprL, -1);

			}
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for (int k = 0; k < data.numScenarios; k++) {
			for (int i = 0; i <data.n; i++) {
				obj.addTerm(data.reward[i]*(1-data.scenarios[k][i]), x[k][i]);
			}
		}

		cplex.addMaximize(obj,"Reward");
		// Optimize model
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Upper Bounder is infeasible!");		
		}
		else{
			// Get solution
			SolutionBilevel sol = new SolutionBilevel(data); // In case we stopped because of time limit
			if(cplex.getStatus() == IloCplex.Status.Optimal || cplex.getStatus() == IloCplex.Status.Feasible) {
				for (int i = 0; i <data.n; i++) {
					sol.probed[i] = (int)Math.round(cplex.getValue(probeVar[i]));
					for (int k = 0; k < data.numScenarios; k++) {
						sol.vl[k][i] = (int)Math.round(cplex.getValue(x[k][i]));
					}
				}
			}
			sol.objl = round(cplex.getBestObjValue());
			sol.objl = sol.objl/(data.numScenarios+0.0);
			sample.add(sol);
			// Update UB	
			UB = Math.min(UB, sol.objl); // In case we stopped because of time limit
			//sol.print();
			//sol.printRewards(data);
		}
		cplex.clearModel();
		cplex=null;
		System.gc();
	}

	public double REHPPScenario(DataHandler data, AlgorithmHandler a, int k, double atime, BBNode node) throws IloException {
		cplex = new IloCplex();
		double timeRemaining = a.timeLimit - (System.currentTimeMillis()-atime)/1000;
		cplex.setParam(IloCplex.Param.TimeLimit, timeRemaining);
		cplex.setParam(IloCplex.IntParam.MIPDisplay, 0);
		cplex.setOut(null);

		IloNumVar[] probeVar = new IloNumVar[data.n];   				// 1 if node is probed
		IloNumVar[] x = new IloNumVar[data.n];   						// 1 if node is visited in scenario k
		IloNumVar[] binaryMultVar = new IloNumVar[data.n];   			// x_i*probe_ki

		//Create variables 
		for (int i = 0; i < data.n; i++) {
			probeVar[i] = cplex.numVar(node.probeLB[i], node.probeUB[i], IloNumVarType.Bool);
		}

		for (int i = 0; i < data.n; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}

		for (int i = 0; i < data.n; i++) {
			binaryMultVar[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}


		// KP constraints
		IloLinearNumExpr expr0 = cplex.linearNumExpr();			
		for (int j = 0; j <data.n; j++) {
			expr0.addTerm(data.weight[j], x[j]);
		}
		cplex.addLe(expr0, data.rhs);

		// Probing budget
		IloLinearNumExpr exprB = cplex.linearNumExpr();
		for (int i = 0; i < data.n; i++) {
			exprB.addTerm(1, probeVar[i]);		
		}
		cplex.addLe(exprB, data.budget);
		//////////////////////////////////////////////////////////////////////////////////////////////

		// Add obj value constraints: Leader Sol Obj >= Sample Sol Obj 
		for (int j = 0; j < sample.size(); j++) {
			SolutionBilevel sol = sample.get(j);
			IloLinearNumExpr expr = cplex.linearNumExpr();
			// Leader Sol obj
			double rhs = 0;
			for (int i = 0; i <data.n; i++) {
				//RHS
				expr.addTerm(data.reward[i]*(1-data.scenarios[k][i]), binaryMultVar[i]);
				expr.addTerm(data.reward[i]*(1-data.probFailure[i]), x[i]);
				expr.addTerm(-data.reward[i]*(1-data.probFailure[i]), binaryMultVar[i]);
				//LHS
				expr.addTerm(-data.reward[i]*(1-data.scenarios[k][i])*sol.xhatf[k][i], probeVar[i]);
				expr.addTerm(data.reward[i]*(1-data.probFailure[i])*sol.xhatf[k][i], probeVar[i]);
				rhs += data.reward[i]*(1-data.probFailure[i])*sol.xhatf[k][i];
			}
			cplex.addGe(expr, rhs-tol);
		}

		//Binary multiplication envelope
		for (int i = 0; i < data.n; i++) {
			IloLinearNumExpr exprA = cplex.linearNumExpr();
			exprA.addTerm(1, binaryMultVar[i]);
			exprA.addTerm(-1, x[i]);
			cplex.addLe(exprA, 0);

			IloLinearNumExpr exprC = cplex.linearNumExpr();
			exprC.addTerm(1, binaryMultVar[i]);
			exprC.addTerm(-1, probeVar[i]);
			cplex.addLe(exprC, 0);

			IloLinearNumExpr exprL = cplex.linearNumExpr();
			exprL.addTerm(1, binaryMultVar[i]);
			exprL.addTerm(-1, probeVar[i]);
			exprL.addTerm(-1, x[i]);
			cplex.addGe(exprL, -1);
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for (int i = 0; i <data.n; i++) {
			obj.addTerm(data.reward[i]*(1-data.scenarios[k][i]), x[i]);
		}

		cplex.addMaximize(obj,"Reward");
		// Optimize model
		cplex.solve();
		double ub = 999999;
		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Upper Bounder is infeasible!");		
		}
		else{
			// Get solution
			SolutionBilevel sol = new SolutionBilevel(data); // In case we stopped because of time limit
			if(cplex.getStatus() == IloCplex.Status.Optimal || cplex.getStatus() == IloCplex.Status.Feasible) {
				for (int i = 0; i <data.n; i++) {
					sol.probed[i] = (int)Math.round(cplex.getValue(probeVar[i]));
					sol.vl[k][i] = (int)Math.round(cplex.getValue(x[i]));
				}
			}
			sol.objl = round(cplex.getBestObjValue()); // In case we stopped because of time limit
			sample.add(sol);
			// Update UB	
			ub = sol.objl;
			//sol.print();
			//sol.printRewards(data);

		}
		cplex.clearModel();
		cplex=null;
		System.gc();
		return ub;
	}

	public void solveFollower(DataHandler data, int nsol, AlgorithmHandler a) throws UnknownObjectException, IloException {
		//Solve for each scenario
		double expectedReward = 0;
		for (int k = 0; k < data.numScenarios; k++) {
			expectedReward += solveSubScenario(data, data.scenarios[k], k, sample.get(nsol), a);
		}
		expectedReward = expectedReward/(data.numScenarios+0.0);

		sample.get(nsol).actualReward = expectedReward;
		sample.get(nsol).estimatedReward = sample.get(nsol).estimatedReward/(data.numScenarios+0.0);

		sample.get(nsol).print();
		//sample.get(nsol).printRewards(data);
	}


	public void sampleHPPCPLEX(DataHandler data, AlgorithmHandler a) throws IloException {
		cplex.setParam(IloCplex.IntParam.MIPDisplay, 0);
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

		cplex.addMaximize(obj,"EstimatedReward");
		cplex.solve();

		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		SAMPLE is infeasible!");
		}
		else {
			// Initial solution with no probing, leader and follower agree! 
			noProbeSol = new SolutionBilevel(data);

			for (int i = 0; i < data.n; i++) {
				if(cplex.getValue(x[i]) > 0.1) { 
					for (int j = 0; j < data.numScenarios; j++) {
						noProbeSol.xhatf[j][i] = 1;						
					}
				}		
			}
			// Estimated and actual reward are the same when no probing
			noProbeSol.estimatedReward = cplex.getObjValue();
			noProbeSol.actualReward = 0;
			for (int i = 0; i < data.n; i++) {
				for (int k = 0; k < data.numScenarios; k++) {
					if(data.scenarios[k][i] == 0 && noProbeSol.xhatf[k][i] >= 0.9) {noProbeSol.actualReward+=data.reward[i];}					
				}
			}

			noProbeSol.actualReward = noProbeSol.actualReward/(data.numScenarios+0.0);
			sample.add(noProbeSol);
			//noProbeSol.print();
			LB = noProbeSol.actualReward; //Update initial LB
		}


		cplex.clearModel();
	}

	public double solveSubScenario(DataHandler data, int[] scenario, int s, SolutionBilevel sol, AlgorithmHandler a) throws UnknownObjectException, IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 3600);
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
			if(sol.probed[i] == 0) {obj.addTerm((1-data.probFailure[i])*data.reward[i], x[i]);}
			else {
				if(scenario[i] == 0) {obj.addTerm(data.reward[i], x[i]);}
				if(scenario[i] == 1) {obj.addTerm(0, x[i]);}
			}
		}

		cplex.addMaximize(obj,"Estimated Reward");
		cplex.solve();


		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Subproblem is infeasible!");
			return -999999;
		}
		else {
			//System.out.println("subproblem estimated reward: "+ (cplex.getObjValue()));
			sol.estimatedReward += cplex.getObjValue();
			double actualReward = 0;
			for (int i = 0; i < data.n; i++) {
				if(scenario[i] == 0 && cplex.getValue(x[i]) >= 0.1) {actualReward+=data.reward[i];}
				if(cplex.getValue(x[i]) >= 0.1) { sol.xhatf[s][i] = 1;}
			}

			// Check if follower is giving an alternative optimal solution to Leader
			double estimatedRL = 0;
			for (int i = 0; i < data.n; i++) {
				if(sol.probed[i] == 1) {estimatedRL += sol.vl[s][i]*data.reward[i]*(1 - scenario[i]);}
				else {estimatedRL += sol.vl[s][i]*data.reward[i]*(1-data.probFailure[i]);}
			}

			if(estimatedRL >= cplex.getObjValue()-tol) { // If solutions are alternative optimal, keep the one from the leader!!! 
				for (int i = 0; i < data.n; i++) {
					sol.xhatf[s][i] = sol.vl[s][i];
				}
				actualReward = actualReward(data, scenario, sol.xhatf[s]);
			}

			cplex.clearModel();
			cplex=null;
			System.gc();
			return actualReward;
		}

	}

	public int[] solveSubScenarioH(DataHandler data, int[] scenario, AlgorithmHandler a) throws UnknownObjectException, IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, 3600);
		cplex.setOut(null);

		IloNumVar[] x = new IloNumVar[data.n];

		//Create variables 
		for (int i = 0; i < data.n; i++) {
			x[i] = cplex.numVar(0, 1, IloNumVarType.Bool);		
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
			if(a.probed[i] == 0) {obj.addTerm((1-data.probFailure[i])*data.reward[i], x[i]);}
			else {
				if(scenario[i] == 0) {obj.addTerm(data.reward[i], x[i]);}
				if(scenario[i] == 1) {obj.addTerm(0, x[i]);}
			}
		}

		cplex.addMaximize(obj,"Estimated Reward");
		cplex.solve();


		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Subproblem HEURISTIC is infeasible!");
			return null;
		}
		else {
			//System.out.println("subproblem HEURISTIC estimated reward: "+ (cplex.getObjValue()));
			int[] xhat = new int[data.n];
			for (int i = 0; i < data.n; i++) {
				if(cplex.getValue(x[i]) >= 0.9) { xhat[i] = 1;}
			}
			cplex.clearModel();
			cplex.end();
			cplex=null;
			System.gc();
			return xhat;

		}

	}

	public void getIncumbentSol(DataHandler data, int[] probeStar) {

		for (int i = 0; i < sample.size(); i++) {
			if(sample.get(i).actualReward>=LB-tol){
				xstar.copy(sample.get(i), data);
				i=sample.size()+1;
			}
		}
		for (int i = 0; i < data.n; i++) {
			probeStar[i] = xstar.probed[i];
		}

	}

	public void getIncumbentSol(DataHandler data, double lb, int q) {

		for (int i = q; i < sample.size(); i++) {
			if(sample.get(i).actualReward>=lb-tol){
				xstar.copy(sample.get(i), data);
				i=sample.size()+1;
			}
		}

	}



	public double[] initializeBB(DataHandler data, double atime, AlgorithmHandler a) throws IloException {

		int[][] probeVar = new int[data.numScenarios][data.n]; // One copy for each scenario
		double[] probeVarMean =  new double[data.n]; // Mean of the different solutions
		double rootBound = 0;
		int[] lb = new int[data.n];
		int[] ub = new int[data.n];
		for (int j = 0; j < data.n; j++) {
			lb[j] = 0;
			ub[j] = 1;
		}
		BBNode node = new BBNode ( rootBound, lb, ub, 0); // Root node
		boolean timeLimitStop = false;

		for (int k = 0; k < data.numScenarios; k++) { // Solve one bilevel problem per scenario
			if((System.currentTimeMillis()-atime)/1000<a.timeLimit) {	
				SolutionBilevel localSol = solveBilevelScenario(data, k, a, atime, node);
				probeVar[k] = localSol.probed;
				for (int i = 0; i < data.n; i++) {
					probeVarMean[i]+=probeVar[k][i];
				}
				rootBound += localSol.actualReward;
				//System.out.println("********************Probed Local: "+Arrays.toString(probeVar[k])+" local obj: "+localSol.actualReward );
			}
			else {
				timeLimitStop = true;
			}
		}
		if(!timeLimitStop) {
			rootBound = rootBound/(data.numScenarios+0.0);
			System.out.println("ROOT BOUND: "+rootBound);
			UB = rootBound;
			node.bound = rootBound;

			// Create first two branches
			double mostFractional = 0;
			int branchIndex = -1;
			for (int i = 0; i < data.n; i++) {
				probeVarMean[i] = probeVarMean[i]/(data.numScenarios+0.0);
				if(Math.min(probeVarMean[i], 1-probeVarMean[i]) > mostFractional) {
					mostFractional = Math.min(probeVarMean[i], 1-probeVarMean[i]);
					branchIndex = i;
				}
			}
			System.out.println("WEPA PROBE MEAN: "+Arrays.toString(probeVarMean));
			if(mostFractional <= tol) {//Root node is integer, we are done!!!!
				System.out.println("Root node is INTEGER!!!!");
				LB = rootBound;
			}
			else {
				System.out.println("Most fractional value: "+mostFractional+" branch var index: "+branchIndex);
				tree = new ArrayList<BBNode>();
				BBNode leftn = new BBNode(node.bound, node.probeLB, node.probeUB, node.numBranches);
				BBNode rightn = new BBNode(node.bound, node.probeLB, node.probeUB, node.numBranches);

				leftn.probeUB[branchIndex] = 0; 
				rightn.probeLB[branchIndex] = 1; 
				tree.add(rightn);
				tree.add(leftn);				
			}

		}
		return probeVarMean;
	}



	private SolutionBilevel solveBilevelScenario(DataHandler data, int k, AlgorithmHandler a, double atime, BBNode node) throws IloException {
		//sample = new ArrayList<SolutionBilevel>();
		int startingSampleSize = sample.size(); 
		int nsol = sample.size();
		double lb = 0;
		double ub = 999999;

		//Start the main loop
		while(ub - lb > tol && (System.currentTimeMillis()-atime)/1000<a.timeLimit){	
			// Get sample size
			nsol = sample.size();
			//Solve REHPP
			ub = REHPPScenario(data, a, k, atime, node);
			// Solve follower's subproblem
			if((System.currentTimeMillis()-atime)/1000<a.timeLimit) {
				//Solve for each scenario
				double expectedReward = solveSubScenario(data, data.scenarios[k], k, sample.get(nsol), a);
				sample.get(nsol).actualReward = expectedReward;
				//sample.get(nsol).print();
				//sample.get(nsol).printRewards(data);

			}
			// Update LB
			if(sample.get(nsol).actualReward> lb){ 
				lb = sample.get(nsol).actualReward; 
			};
			//System.out.println("Scenario : "+k+" UB: "+ub+" LB: "+lb+" Sample size: "+nsol);
		}	

		// Print the results
		SolutionBilevel localSol = new SolutionBilevel(data);
		getIncumbentSol(data, lb, startingSampleSize);
		localSol.copy(xstar, data);
		return localSol;


	}



	public void runBB(DataHandler data, double atime, AlgorithmHandler a) throws IloException {
		double timeUsed = (System.currentTimeMillis()-atime)/1000;
		int iter = 1;
		while(tree.size() > 0 && (UB-LB)/(Math.abs(UB)) > 0.001 && timeUsed < a.timeLimit) {
			BBNode n = tree.get(0);
			tree.remove(0);
			numNodesExplored++;
			solveNode(data, atime, a ,n);

			timeUsed = (System.currentTimeMillis()-atime)/1000;
			if(tree.size()>0 && timeUsed < a.timeLimit) {
				UB = Math.min(UB, tree.get(0).bound);
				if(tree.get(0).bound < LB) { // We are done!!!!
					UB = LB;
				}
			}
			else if (tree.size() == 0 && timeUsed < a.timeLimit) { // We are done!!! 
				UB = LB;
			}
			iter++;
			timeUsed = (System.currentTimeMillis()-atime)/1000;

			if(iter%1==0) {
				System.out.println("NODES EXPLORED "+iter+" LB = "+LB+" UB = "+UB+" GAP: "+round( (UB-LB)/(Math.abs(UB)) ) );
				System.out.println("TIME USED: "+timeUsed);
			}

		}


	}



	private void solveNode(DataHandler data, double atime, AlgorithmHandler a, BBNode node) throws IloException {

		//System.out.println("SOLVING NODE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		//node.print();
		int[][] probeVar = new int[data.numScenarios][data.n]; // One copy for each scenario
		double[] probeVarMean =  new double[data.n]; // Mean of the different solutions
		boolean timeLimitStop = false;
		double localBound = 0;

		for (int k = 0; k < data.numScenarios; k++) { // Solve one bilevel problem per scenario
			if((System.currentTimeMillis()-atime)/1000<a.timeLimit) {	
				SolutionBilevel localSol = solveBilevelScenario(data, k, a, atime, node);
				probeVar[k] = localSol.probed;
				for (int i = 0; i < data.n; i++) {
					probeVarMean[i]+=probeVar[k][i];
				}
				localBound += localSol.actualReward;
				//System.out.println("********************Probed Local: "+Arrays.toString(probeVar[k])+" local obj: "+localSol.actualReward );
			}
			else {
				timeLimitStop = true;
			}
		}
		if(!timeLimitStop) {
			localBound = localBound/(data.numScenarios+0.0);
			//System.out.println("LOCAL BOUND: "+localBound);
			node.bound = localBound;

			// Create next two branches
			double mostFractional = 0;
			int branchIndex = -1;
			for (int i = 0; i < data.n; i++) {
				probeVarMean[i] = probeVarMean[i]/(data.numScenarios+0.0);
				if(Math.min(probeVarMean[i], 1-probeVarMean[i]) > mostFractional) {
					mostFractional = Math.min(probeVarMean[i], 1-probeVarMean[i]);
					branchIndex = i;
				}
			}
			//System.out.println("WEPA PROBE MEAN: "+Arrays.toString(probeVarMean));
			if(mostFractional <= tol) {//Root node is integer, we are done!!!!
				System.out.println("Local Solution is INTEGER!!!!");
				LB = Math.max(LB, localBound);
			}
			else if(node.bound > LB){
				//System.out.println("Most fractional value: "+mostFractional+" branch var index: "+branchIndex);
				BBNode leftn = new BBNode(node.bound, node.probeLB, node.probeUB, node.numBranches);
				BBNode rightn = new BBNode(node.bound, node.probeLB, node.probeUB, node.numBranches);

				leftn.probeUB[branchIndex] = 0; 
				rightn.probeLB[branchIndex] = 1; 
				tree.add(rightn);
				tree.add(leftn);


				Collections.sort(tree, Collections.reverseOrder());
				//printTree();

			}

		}




	}



	private void printTree() {
		System.out.println("*********************************************************************************************************");
		System.out.println("BB TREE SIZE "+tree.size()+" BOUND: "+tree.get(0).bound+" Branches: "+tree.get(0).numBranches);
		/*for (int i = 0; i < tree.size(); i++) {
			System.out.println("BB TREE INDEX "+i+" BOUND: "+tree.get(i).bound+" Branches: "+tree.get(i).numBranches);
		}*/
		System.out.println("*********************************************************************************************************");
	}

	// Heuristic used in the BB tree
	public void primalHeuristic(DataHandler data, double atime, AlgorithmHandler a, double[] meanProbe) throws UnknownObjectException, IloException {
		// Sort the mean probe vector and select the top B
		List<Item> order=new ArrayList<>();
		//System.out.println("Primal heuristic mean probe "+Arrays.toString(meanProbe));

		for (int i = 0; i < data.n; i++) {
			order.add(new Item(i, meanProbe[i] ));
		}
		Collections.sort(order);

		a.probed = new int[data.n];
		for (int i = 0; i < data.budget; i++) {
			//System.out.println("Selected probe: "+order.get(i).pos);
			a.probed[order.get(i).pos] = 1; 
		}
		System.out.println("Primal heuristic probe "+Arrays.toString(a.probed));

		// Fix the probe decisions and solve for each scenario
		solveForFixedProbe(data, atime, a);

	}

	// Heuristic sorting reward
	public void greedyHeuristic(DataHandler data, double atime, AlgorithmHandler a) throws UnknownObjectException, IloException {
		// Sort the reward vector and select the top B
		List<Item> order=new ArrayList<>();

		for (int i = 0; i < data.n; i++) {
			order.add(new Item(i, data.reward[i] ));
		}
		Collections.sort(order);

		double bestObj = -999999;
		HashMap<String, Double> probeMap = new HashMap<String, Double> ();	//Maps the probes we already visited
		for (int q = 0; q < 10; q++) { //Main loop
			a.probed = new int[data.n];
			int numProbed = 0;
			for (int i = 0; i < data.n; i++) {
				if(data.r.nextDouble() <= 0.75) {
					a.probed[order.get(i).pos] = 1;
					numProbed++;
					//System.out.println("PROBED POSITION "+order.get(i).pos+" reward: "+data.reward[order.get(i).pos]);
				}
				if(numProbed == data.budget) {i = data.n;} // No more budget!!!
			}
			//System.out.println("Greedy heuristic reward probe "+Arrays.toString(a.probed));
			if(probeMap.containsKey(Arrays.toString(a.probed))) {
				//System.out.println("Repeated probe!!!");
				q--;
			}
			else {
				// Fix the probe decisions and solve for each scenario
				double objH = solveForFixedProbe(data, atime, a);
				probeMap.put(Arrays.toString(a.probed), objH);
				//System.out.println("OBJ: "+objH);
				if(objH > bestObj) {
					bestObj = objH;
				}
			}
		}
		System.out.println("Greedy heuristic reward probe "+Arrays.toString(a.probed));
		data.greedyObj1 = bestObj;
		System.out.println("Heuristic Obj: "+data.greedyObj1);
		System.out.println();

	}

	// Heuristic sorting estimated reward
	public void greedyHeuristicPI(DataHandler data, double atime, AlgorithmHandler a) throws UnknownObjectException, IloException {
		
		// Sort the reward vector according to components used in perfect info
		int[] setQ = findSetQ(data, atime, a); //See section 5 of the paper!
		int countNonZero = 0;
		List<Item> order=new ArrayList<>();
		for (int i = 0; i < data.n; i++) {
			order.add(new Item(i, setQ[i]));
			if(setQ[i] >= 1) {
				countNonZero++;
			}
		}
		Collections.sort(order);

		
		System.out.println("MAXIMUM BUDGET NEEDED FOR PI PEFORMANCE: "+countNonZero);
		data.budgetUB = countNonZero;

		double bestObj = -999999;
		HashMap<String, Double> probeMap = new HashMap<String, Double> ();	//Maps the probes we already visited
		for (int q = 0; q < 10; q++) { //Main loop
			a.probed = new int[data.n];
			int numProbed = 0;
			for (int i = 0; i < data.n; i++) {
				if(data.r.nextDouble() <= 0.75) {
					a.probed[order.get(i).pos] = 1;
					numProbed++;
					//System.out.println("PROBED POSITION "+order.get(i).pos+" count: "+count[order.get(i).pos]);
				}
				if(numProbed == data.budget) {i = data.n;} // No more budget!!!
			}
			//System.out.println("Greedy heuristic expected reward probe "+Arrays.toString(a.probed));
			if(probeMap.containsKey(Arrays.toString(a.probed))) {
				//System.out.println("Repeated probe!!!");
				q--;
			}
			else {
				// Fix the probe decisions and solve for each scenario
				double objH = solveForFixedProbe(data, atime, a);
				probeMap.put(Arrays.toString(a.probed), objH);
				//System.out.println("OBJ: "+objH);
				if(objH > bestObj) {
					bestObj = objH;
				}
			}
		}

		data.greedyObjPI = bestObj;
		System.out.println("Heuristic Obj: "+data.greedyObjPI);
		System.out.println();


	}
	
	private int[] findSetQ(DataHandler data, double atime, AlgorithmHandler a) throws UnknownObjectException, IloException {
		// Set the initial Q
		int[] setQ = new int[data.n];
		a.probed = new int[data.n];
		solvePerfectInfo(data, a);
		for (int i = 0; i < data.n; i++) {
			for (int k = 0; k < data.numScenarios; k++) {
				if(perfectInfoSol.xhatf[k][i] >= 0.1) {
					setQ[i]++;
					a.probed[i] = 1;
				}
			}
		}
		
		boolean stop = false;
		while(!stop) {
			stop = true;
			// Solve the follower problem for each scenario and check if the estimated rewards match
			for (int k = 0; k < data.numScenarios; k++) {
					int[] xhat = solveSubScenarioH(data, data.scenarios[k], a);
					double eCost = estimatedReward(data, data.scenarios[k], a.probed, xhat);
					double eCostPI = estimatedReward(data, data.scenarios[k], a.probed, perfectInfoSol.xhatf[k]);
					
					if(Math.abs(eCostPI-eCost) > 0.0001) {
						System.out.println("Estimated cost doesn't match: "+eCost+" vs "+eCostPI);
						stop = false;
						// Find an index not in Q where the solutions won't agree and add it to set Q
						for (int j = 0; j < data.n; j++) {
							if(setQ[j] == 0 && (xhat[j] != perfectInfoSol.xhatf[k][j])) {
								System.out.println("ADDING INDEX "+j+" TO SET Q");
								setQ[j] = 1;
								a.probed[j] = 1;
								j = data.n;
							}
						}
						k = data.numScenarios; //Stop outer for loop
					}
					
			}
			
		}
		
		
		return setQ;
	}

	

	// Heuristic: find probe that maximizes estimated reward
	public void heuristicER(DataHandler data, double atime, AlgorithmHandler a) throws UnknownObjectException, IloException {
		double bestObj = -999999;
		for (int k = 0; k < data.numScenarios; k++) {
			a.probed = new int[data.n];
			a.probed = findBestProbeERHeuristic(data, k, a);
			System.out.println("Heuristic ER probe "+Arrays.toString(a.probed));
			// Fix the probe decisions and solve for each scenario
			double objH = solveForFixedProbe(data, atime, a);
			System.out.println("OBJ H: "+objH);
			if(objH > bestObj) {
				bestObj = objH;
			}

		}
		data.heuristicER = bestObj;
		System.out.println("Heuristic Obj: "+data.heuristicER);
		System.out.println();

	}

	public void heuristicER2(DataHandler data, double atime, AlgorithmHandler a) throws UnknownObjectException, IloException {
		int[][] probed = new int[data.n][20];
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, a.timeLimit/2);
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 0.05);
		//cplex.setParam(IloCplex.IntParam.MIPDisplay, 0);
		//cplex.setOut(null);

		IloNumVar[][] x = new IloNumVar[data.numScenarios][data.n];   	
		IloNumVar[] probeVar = new IloNumVar[data.n];   										// 1 if node is probed
		IloNumVar[][] binaryMultVar = new IloNumVar[data.numScenarios][data.n];   				// x_ki*probe_ki

		//Create variables 
		for (int i = 0; i < data.n; i++) {
			for (int k = 0; k < data.numScenarios; k++) {
				x[k][i]= cplex.numVar(0, 1, IloNumVarType.Bool);		
			}
		}
		for (int i = 0; i < data.n; i++) {
			probeVar[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		for (int i = 0; i < data.n; i++) {
			for (int k = 0; k < data.numScenarios; k++) {
				binaryMultVar[k][i] = cplex.numVar(0, 1, IloNumVarType.Bool);
			}
		}

		// KP constraints
		for (int k = 0; k < data.numScenarios; k++) {
			IloLinearNumExpr expr0 = cplex.linearNumExpr();			
			for (int j = 0; j <data.n; j++) {
				expr0.addTerm(data.weight[j], x[k][j]);
			}
			cplex.addLe(expr0, data.rhs);
		}

		// Probing budget
		IloLinearNumExpr exprB = cplex.linearNumExpr();
		for (int i = 0; i < data.n; i++) {
			exprB.addTerm(1, probeVar[i]);		
		}
		cplex.addLe(exprB, data.budget);


		//////////////////////////////////////////////////////////////////////////////////////////////
		//Binary multiplication envelope
		for (int i = 0; i < data.n; i++) {
			for (int k = 0; k < data.numScenarios; k++) {
				IloLinearNumExpr exprA = cplex.linearNumExpr();
				exprA.addTerm(1, binaryMultVar[k][i]);
				exprA.addTerm(-1, x[k][i]);
				cplex.addLe(exprA, 0);

				IloLinearNumExpr exprC = cplex.linearNumExpr();
				exprC.addTerm(1, binaryMultVar[k][i]);
				exprC.addTerm(-1, probeVar[i]);
				cplex.addLe(exprC, 0);

				IloLinearNumExpr exprL = cplex.linearNumExpr();
				exprL.addTerm(1, binaryMultVar[k][i]);
				exprL.addTerm(-1, probeVar[i]);
				exprL.addTerm(-1, x[k][i]);
				cplex.addGe(exprL, -1);

			}
		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for (int k = 0; k < data.numScenarios; k++) {
			for (int i = 0; i <data.n; i++) {
				obj.addTerm(data.reward[i]*(1-data.scenarios[k][i]), binaryMultVar[k][i]);
				obj.addTerm(data.reward[i]*(1-data.probFailure[i]), x[k][i]);
				obj.addTerm(-data.reward[i]*(1-data.probFailure[i]), binaryMultVar[k][i]);
			}
		}

		cplex.addMaximize(obj,"estimatedReward");
		// Optimize model
		for (int q = 0; q < 10; q++) {
			cplex.solve();
			if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
				System.out.println("		Heuristic is infeasible!");		
			}
			else{
				// Get solution
				if(cplex.getStatus() == IloCplex.Status.Optimal || cplex.getStatus() == IloCplex.Status.Feasible) {
					for (int i = 0; i <data.n; i++) {
						probed[i][q] = (int)Math.round(cplex.getValue(probeVar[i]));
					}
					// Cut solution (at least 1 different probes!)
					int count = 0;
					IloLinearNumExpr exprC = cplex.linearNumExpr();
					for (int i = 0; i < data.n; i++) {
						if(probed[i][q] == 1) {exprC.addTerm(1, probeVar[i]); count++;}		
					}
					cplex.addLe(exprC, count-1);
					
				}
			}
			
		}
		cplex.clearModel();
		cplex=null;
		System.gc();

		double bestObj = -999999;
		for (int q = 0; q < 10; q++) {
			a.probed = new int[data.n];
			for (int i = 0; i < data.n; i++) {
				a.probed[i] = probed[i][q];
			}
			System.out.println("Heuristic ER2 probe "+Arrays.toString(a.probed));
			// Fix the probe decisions and solve for each scenario
			double objH = solveForFixedProbe(data, atime, a);
			System.out.println("OBJ H: "+objH);
			if(objH > bestObj) {
				bestObj = objH;
			}

		}
		data.heuristicER = bestObj;
		System.out.println("Heuristic Obj: "+data.heuristicER);
		System.out.println();

	}

	private int[] findBestProbeERHeuristic(DataHandler data, int k, AlgorithmHandler a) throws IloException {
		cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.TimeLimit, a.timeLimit/data.numScenarios);
		cplex.setParam(IloCplex.IntParam.MIPDisplay, 0);
		cplex.setOut(null);

		IloNumVar[] x = new IloNumVar[data.n];   	
		IloNumVar[] probeVar = new IloNumVar[data.n];   				// 1 if node is probed
		IloNumVar[] binaryMultVar = new IloNumVar[data.n];   	// x_ki*probe_ki

		//Create variables 
		for (int i = 0; i < data.n; i++) {
			x[i]= cplex.numVar(0, 1, IloNumVarType.Bool);		
		}
		for (int i = 0; i < data.n; i++) {
			probeVar[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}
		for (int i = 0; i < data.n; i++) {
			binaryMultVar[i] = cplex.numVar(0, 1, IloNumVarType.Bool);
		}

		// KP constraints
		IloLinearNumExpr expr0 = cplex.linearNumExpr();			
		for (int j = 0; j <data.n; j++) {
			expr0.addTerm(data.weight[j], x[j]);
		}
		cplex.addLe(expr0, data.rhs);

		// Probing budget
		IloLinearNumExpr exprB = cplex.linearNumExpr();
		for (int i = 0; i < data.n; i++) {
			exprB.addTerm(1, probeVar[i]);		
		}
		cplex.addLe(exprB, data.budget);


		//////////////////////////////////////////////////////////////////////////////////////////////
		//Binary multiplication envelope
		for (int i = 0; i < data.n; i++) {
			IloLinearNumExpr exprA = cplex.linearNumExpr();
			exprA.addTerm(1, binaryMultVar[i]);
			exprA.addTerm(-1, x[i]);
			cplex.addLe(exprA, 0);

			IloLinearNumExpr exprC = cplex.linearNumExpr();
			exprC.addTerm(1, binaryMultVar[i]);
			exprC.addTerm(-1, probeVar[i]);
			cplex.addLe(exprC, 0);

			IloLinearNumExpr exprL = cplex.linearNumExpr();
			exprL.addTerm(1, binaryMultVar[i]);
			exprL.addTerm(-1, probeVar[i]);
			exprL.addTerm(-1, x[i]);
			cplex.addGe(exprL, -1);

		}

		// ADD OBJECTIVE FUNCTION	
		IloLinearNumExpr obj = cplex.linearNumExpr();
		for (int i = 0; i <data.n; i++) {
			obj.addTerm(data.reward[i]*(1-data.scenarios[k][i]), binaryMultVar[i]);
			obj.addTerm(data.reward[i]*(1-data.probFailure[i]), x[i]);
			obj.addTerm(-data.reward[i]*(1-data.probFailure[i]), binaryMultVar[i]);
		}
		cplex.addMaximize(obj,"EstimatedReward");
		// Optimize model
		cplex.solve();
		int[] probed = new int[data.n];
		if(cplex.getStatus() == IloCplex.Status.Infeasible || cplex.getStatus() == IloCplex.Status.InfeasibleOrUnbounded){
			System.out.println("		Heuristic ER is infeasible!");		
		}
		else{
			// Get solution

			if(cplex.getStatus() == IloCplex.Status.Optimal || cplex.getStatus() == IloCplex.Status.Feasible) {
				for (int i = 0; i <data.n; i++) {
					probed[i] = (int)Math.round(cplex.getValue(probeVar[i]));
				}
			}
			System.out.println("Heuristic best estimated reward: "+round(cplex.getBestObjValue()));
		}
		cplex.clearModel();
		cplex=null;
		System.gc();
		return probed;

	}



	private double solveForFixedProbe(DataHandler data, double atime, AlgorithmHandler a) throws UnknownObjectException, IloException {
		double expectedReward = 0;
		HashMap<String, int[]> probeResults = new HashMap<String, int[]> (); 
		for (int k = 0; k < data.numScenarios; k++) {
			// Store current info from the probe
			int [] pResult = new int[data.n]; 
			for (int i = 0; i < data.n; i++) {// 0 not probed, -1 no jammer confirmed, 1 jammer confirmed
				if(a.probed[i] == 1) {
					if(data.scenarios[k][i] == 0) {pResult[i] = -1;}
					if(data.scenarios[k][i] == 1) {pResult[i] = 1;}
				}
			}
			if(probeResults.containsKey(Arrays.toString(pResult))) {
				int[] xhat = probeResults.get(Arrays.toString(pResult));
				expectedReward += actualReward(data, data.scenarios[k], xhat);
			}

			else {
				int[] xhat = solveSubScenarioH(data, data.scenarios[k], a);
				expectedReward += actualReward(data, data.scenarios[k], xhat);
				probeResults.put(Arrays.toString(pResult), xhat);
			}

		}
		expectedReward = expectedReward/(data.numScenarios+0.0);

		// Update best solution
		if(expectedReward > LB) {
			LB = expectedReward;
			for (int i = 0; i <data.n; i++) {
				a.probeStar[i] = a.probed[i];
			}
			System.out.println("*************************************************NEW BEST SOLUTION EXPECTED REWARD: "+LB);
		}

		//System.out.println("EXPECTED ACTUAL REWARD: "+expectedReward);
		return expectedReward;

	}

	private double actualReward(DataHandler data, int[] scenario, int[] x) {
		double actualReward = 0;
		for (int i = 0; i < data.n; i++) {
			if(scenario[i] == 0 && x[i] >= 0.9) {
				actualReward+=data.reward[i]; 
			}
		}

		//System.out.println("subproblem ACTUAL reward: "+actualReward);

		return actualReward;
	}
	
	private double estimatedReward(DataHandler data, int[] scenario, int[] probed, int[] x) {
		double eCost = 0;
		for (int i = 0; i < data.n; i++) {
			if(probed[i] == 0) {eCost += ((1-data.probFailure[i])*data.reward[i]+data.probFailure[i]*0)*x[i];}
			else {
				if(scenario[i] == 0) {eCost += data.reward[i]*x[i];}
				if(scenario[i] == 1) {eCost += 0*x[i];}
			}
		}
		return eCost;
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
			return -Double.compare(val, o.val);
		}

	}



}
