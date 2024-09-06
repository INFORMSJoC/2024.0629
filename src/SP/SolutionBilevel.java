package SP;

import java.util.Arrays;

public class SolutionBilevel {

	// probe variables
	int[] probed;
	// leader's visited variables
	int[][] xhatl;
	// follower's visited variables
	int[][] xhatf;
	// Objective value for the HPP
	double objl;
	// Objective value for the Follower
	double estimatedCost;
	// Objective value
	double actualCost;

	
	SolutionBilevel(DataHandler data){
		probed = new int[data.numArcs];
		xhatl = new int[data.numScenarios][data.numArcs];
		xhatf = new int[data.numScenarios][data.numArcs];
		actualCost = 0;
		objl = 0;
		estimatedCost = 0;
	}

	public void print() {
		/*System.out.println("Leader Obj: "+objl);
		System.out.println("Estimated Reward: "+estimatedReward);
		System.out.println("Actual Reward: "+actualReward);
		System.out.println("probed: "+Arrays.toString(probed));
		for (int k = 0; k < 3; k++) {
			System.out.println("vl: "+Arrays.toString(vl[k]));	
		}
		System.out.println();
		for (int k = 0; k < 3; k++) {
			System.out.println("vf: "+Arrays.toString(vf[k]));	
		}*/
		System.out.println("Probed: ");
		for (int k = 0; k < probed.length; k++) {
			if(probed[k] > 0.1) {System.out.print("z_"+k+" ");}	
		}
		System.out.println();
		
	}
	public void copy(SolutionBilevel seed, DataHandler data) {
		actualCost = seed.actualCost;
		estimatedCost = seed.estimatedCost;
		objl = seed.objl;
		
		for (int i = 0; i < probed.length; i++) {
			probed[i] = seed.probed[i];
		}
		for (int k = 0; k < data.numScenarios; k++) {
			for (int i = 0; i < data.numArcs; i++) {
				xhatf[k][i] = seed.xhatf[k][i];				
			}
		}
		
		
	}



}
