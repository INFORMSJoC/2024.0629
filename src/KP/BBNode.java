package KP;

public class BBNode implements Comparable<BBNode> {

	int[] probeLB;
	int[] probeUB;	
	double bound; 		
	int numBranches;
		
	public BBNode(double _bound, int[] _xlb, int[] _xub, int numB) {
		bound = _bound;
		probeLB = new int[_xlb.length];
		probeUB = new int[_xub.length];
		for (int j = 0; j < _xub.length; j++) {
			probeLB[j] = _xlb[j];
			probeUB[j] = _xub[j];
		}
		numBranches = numB+1;
	}

	public void print() {
		System.out.println("BB NODE BOUND: "+bound);
		for (int i = 0; i < probeLB.length; i++) {
			System.out.print(" x_"+i+" ["+probeLB[i]+" , "+probeUB[i]+"] ");
		}
		System.out.println();
	}
	
	@Override
	public int compareTo(BBNode n) {
		if(this.bound > n.bound) {return 1;}
		else if(this.bound < n.bound) {return -1;}
		else{return 0;}
	}

}
