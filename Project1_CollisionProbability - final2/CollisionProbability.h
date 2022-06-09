#pragma once

#include "define.h"
#include "Matrix.h"

class CollisionProbability
{

private:
	// variable in input file
	int nr;
	int ngauss;
	
	real albedo;
	real jext;

	real *rad;
	real *sigt;
	real *sigs;

	real *sigr;

	real *qi;

	// variable for gauss-jacobi quadrature
	real *x_i;
	real *weight;

	// variable for middle of calculation
	real *delta;	// delta_R_[k]
	

	real **Sij;		// (nr+1)*(nr+1)
	real ***Sijk;		// (nr+1)*(nr+1)*nr		////여기서 k만 index가 0부터 시작된다.
	real **Pij;		// nr*nr
	
	real *vol;		// 1*nr
	real *ci;		// 1*nr
	real **rho;		// (nr+1)*ngauss   

	real gamma;

	// variable for matrix calculation
	real ** A;// (nr);
	real ** Xik;//{ nr };
	real ** Xika;//(nr);

	real * b_y;//(nr, 1);
	real * b_x;//(nr, 1);

	real * Yi;//(nr, 1);
	real * Xiktmp;//(nr, 1);

	real * Ya;//(nr, 1);
	
	real * xk;//(nr, 1);

	real * sourceFlux;//(nr, 1);
	real * J_Flux;//(nr, 1);
	real * Flux;//(nr, 1);

	//
	real *removal;
public:
	string inName;



public:
	CollisionProbability();
	~CollisionProbability();

	void cpman(string inputName);
	void cpman_GaussianQuad(string inputName);
	void cpman_PerturbRemovalXsec(string inputName);

	void manipulateRemovalXsec(real rate);

	void readInput();
	void allocateVaiable();
	
	real calBickeley(real x_in);
	void gaussJacobiQuad(int inGaussNum, real * x_i, real *weight);
	void gaussianQuad(int inGaussNum, real * x_i, real *weight);
	
	void calCPKernels();	// calculate Pij for i,j=1...nr

	void calCPKernels_GaussianQuad();
	
	void solveSystem();


	void luFactorization(real **matA);
	void solveLU(real **matA,real *vec_x, real *vec_b);
	

};

