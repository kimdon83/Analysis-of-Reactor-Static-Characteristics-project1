//Project1 : Collision Probability Calculation for a Circular Cell

#include "define.h"
#include "CollisionProbability.h"
#include "Matrix.h"

int main(void) {

	string name;
	int calType;

	cout << "Collision Probability Method for Circular Cell" << endl;
	cout << "Select the calculation type" << endl;
	cout << "1. Single calculation " << endl;
	cout << "2. Single calculation with gaussian quadrature" << endl;
	cout << "3. Same calculation N time with perturbed removal Xsec" << endl;
	cin >> calType;

	cout << "Select the input" << endl;
	cout << "ex) input1" << endl;
	cout << "ex) input1-2" << endl;
	cout << "ex) input2" << endl;
	cin >> name;

	int n,i;
	if (calType == 3) {
		cout << "how many calculation do you want?" << endl;
		cin >> n;
	}

	//time start
	clock_t start_point, end_point;
	start_point = clock();
	double caltime;
	
	/*name = "input";
	int i = 2;
	name.append(to_string(i));*/
	if (calType == 1) {
		CollisionProbability cp;
		cp.cpman(name);
		//cp.cpman_PerturbRemovalXsec(name);
	}
	else if (calType == 2) {
			CollisionProbability cp;
			cp.cpman_GaussianQuad(name);
		
	}
	else if (calType == 3) {
		for (i = 0; i < n; i++) {
			
			CollisionProbability cp;
			cp.inName = name;
			//cp.cpman(name);
			cp.cpman_PerturbRemovalXsec(name);
		}
	}
	else {
		cout << "Select the calculation Type in the 1,2,3" << endl;
	}
	
	end_point = clock();
	caltime = (end_point - start_point) / (double)(CLOCKS_PER_SEC);

	cout << caltime << "(sec)" << endl;

	system("pause");

	return 0;
}