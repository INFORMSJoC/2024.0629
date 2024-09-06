package KP;

import java.util.Arrays;

public class SolutionBilevel {

	// probe variables
	int[] probed;
	// leader's visited variables
	int[][] vl;
	// follower's visited variables
	int[][] xhatf;
	// Objective value for the HPP
	double objl;
	// Objective value for the Follower
	double estimatedReward;
	// Objective value
	double actualReward;

	
	SolutionBilevel(DataHandler data){
		probed = new int[data.n];
		vl = new int[data.numScenarios][data.n];
		xhatf = new int[data.numScenarios][data.n];
		actualReward = 0;
		objl = 0;
		estimatedReward = 0;
	}

	public void print() {
		System.out.println("Leader Obj: "+objl);
		System.out.println("Estimated Reward: "+estimatedReward);
		System.out.println("Actual Reward: "+actualReward);
		System.out.println("probed: "+Arrays.toString(probed));
		for (int k = 0; k < 3; k++) {
			System.out.println("vl: "+Arrays.toString(vl[k]));	
		}
		System.out.println();
		for (int k = 0; k < 3; k++) {
			System.out.println("vf: "+Arrays.toString(xhatf[k]));	
		}
		
		
	}
	public void copy(SolutionBilevel seed, DataHandler data) {
		actualReward = seed.actualReward;
		estimatedReward = seed.estimatedReward;
		objl = seed.objl;
		
		for (int i = 0; i < probed.length; i++) {
			probed[i] = seed.probed[i];
		}
		for (int k = 0; k < data.numScenarios; k++) {
			for (int i = 0; i < data.n; i++) {
				xhatf[k][i] = seed.xhatf[k][i];				
			}
		}
		
		
	}

	public void printRewards(DataHandler data) {
		for (int k = 0; k < 3; k++) {
			double estimatedRL = 0;
			for (int i = 0; i < data.n; i++) {
				if(probed[i] == 1) {estimatedRL += vl[k][i]*data.reward[i]*(1-data.scenarios[k][i]);}
				else {estimatedRL += vl[k][i]*data.reward[i]*(1-data.probFailure[i]);}
			}
			System.out.println("Scenario "+k+" Leader Estimated Reward: "+estimatedRL);
		}
		System.out.println();
		for (int k = 0; k < 3; k++) {
			double estimatedRF = 0;
			for (int i = 0; i < data.n; i++) {
				if(probed[i] == 1) {estimatedRF += xhatf[k][i]*data.reward[i]*(1-data.scenarios[k][i]);}
				else {estimatedRF += xhatf[k][i]*data.reward[i]*(1-data.probFailure[i]);}
			}
			System.out.println("Scenario "+k+" Follower Estimated Reward: "+estimatedRF);
		}
		
		
	}

}
