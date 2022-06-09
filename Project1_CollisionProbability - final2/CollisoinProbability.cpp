#include "CollisionProbability.h"
#include "Matrix.h"

#include "define.h"

CollisionProbability::CollisionProbability()
{
}
CollisionProbability::~CollisionProbability()
{
	int i, j;

	delete[] rad;
	delete[] sigt;
	delete[] sigs;

	delete[] sigr;
	
	delete[] qi;

	// variable for gauss-jacobi quadrature
	delete[] x_i;
	delete[] weight;

	// variable for middle of calculation
	delete[] delta;	// delta_R_[k]

	for (i = 0; i < nr + 1; i++) {
		for (j = 0; j < nr + 1; j++) {
			delete[] Sijk[i][j];
		}
		delete[] Sijk[i];
	}
	delete[] Sijk;

	for (i = 0; i < nr + 1; i++) {
		delete[] Sij[i];
	}
	delete[] Sij;

	for (i = 0; i < nr; i++) {
		delete[] Pij[i];
	}
	delete[] Pij;

	for (i = 0; i < nr + 1; i++) {
		delete[] rho[i];
	}
	delete[] rho;


	// variable for matrix calculation
	for (i = 0; i < nr ; i++) {
		delete[] A[i];
	}
	delete[] A;

	for (i = 0; i < nr; i++) {
		delete[] Xik[i];
	}
	delete[] Xik;

	for (i = 0; i < nr; i++) {
		delete[] Xika[i];
	}
	delete[] Xika;

	
	delete[] b_y;//(nr, 1);
	delete[] b_x;//(nr, 1);

	delete[] Yi;//(nr, 1);
	delete[] Xiktmp;//(nr, 1);

	delete[] Ya;//(nr, 1);

	delete[] xk;//(nr, 1);

	delete[] sourceFlux;//(nr, 1);
	delete[] J_Flux;//(nr, 1);
	delete[] Flux;//(nr, 1);

	//
	delete[] removal;

}
void CollisionProbability::cpman(string inputName)	//cp로 문제 하나를 푸는 전체과정을 총괄하는 함수
{
	inName = inputName;
	readInput();
	gaussJacobiQuad(ngauss, x_i, weight);
	calCPKernels();

	solveSystem();

}void CollisionProbability::cpman_GaussianQuad(string inputName)	//cp로 문제 하나를 푸는 전체과정을 총괄하는 함수
{
	inName = inputName;
	readInput();
	gaussianQuad(ngauss, x_i, weight);
	calCPKernels_GaussianQuad();

	solveSystem();

}

void CollisionProbability::cpman_PerturbRemovalXsec(string inputName)
{
	{
		inName = inputName;
		readInput();
		gaussJacobiQuad(ngauss, x_i, weight);
		manipulateRemovalXsec(0.03);
		calCPKernels();
		solveSystem();
	}
}

void CollisionProbability::manipulateRemovalXsec(real rate)
{
	srand(time(NULL));

	int i;
	real rn;
	for (i = 0; i < nr; i++) {
		rn = (-1.0+2.*(rand()/RAND_MAX))*rate;
		sigt[i] = (sigt[i] - sigs[i])*(1 + rn) + sigs[i];
	}

}

void CollisionProbability::readInput()
{
	//입력파일이름은 Col..Proba.. 클래스에 string inName 으로 저장되어있다.
	string inputName;
	inputName = inName;
	inputName.append(".txt");

	ifstream fin(inputName);
	string str;

	int i; // for counter

	fin >> str;	// input_set_# 와 빈줄 한줄을 넘어간다. 
	fin >> str;
	nr = atof(str.c_str());
	/*c_str is needed to convert string to const char* previously (the function requires it)*/

	fin >> str;
	ngauss = atof(str.c_str());
	fin >> str;	// !nr...을 넘어간다.

	allocateVaiable();
	// 클래스 자체에서 메모리 주소만 정의되어있던 포인터 변수들을 nr, 즉 크기에 맞게 동적할당해준다.
	// 소멸자에서 반드시 이에 대응하는 delete(deallocate 함수가 있어야 한다.)

	for (i = 0; i < nr; i++) {
		fin >> str;
		rad[i] = atof(str.c_str());
	}
	fin >> str;	// !radius 를 넘어간다.

	for (i = 0; i < nr; i++) {
		fin >> str;
		sigt[i] = atof(str.c_str());
	}
	fin >> str;	// !sigt 를 넘어간다.

	for (i = 0; i < nr; i++) {
		fin >> str;
		sigs[i] = atof(str.c_str());
	}
	fin >> str;	// !sigs 를 넘어간다.

	for (i = 0; i < nr; i++) {
		fin >> str;
		qi[i] = atof(str.c_str());
	}
	fin >> str;	// !qi 를 넘어간다.

	fin >> str;
	albedo = atof(str.c_str());
	fin >> str;
	jext = atof(str.c_str());

	fin.close();
}

void CollisionProbability::allocateVaiable()
{
	int i, j;
	rad = new real[nr];
	memset(rad, 0, sizeof(real)*(nr));
	sigt = new real[nr];
	memset(sigt, 0, sizeof(real)*(nr));
	sigs = new real[nr];
	memset(sigs, 0, sizeof(real)*(nr));

	sigr = new real[nr];
	memset(sigs, 0, sizeof(real)*(nr));


	qi = new real[nr];
	memset(qi, 0, sizeof(real)*(nr));
	x_i = new real[nr];
	memset(x_i, 0, sizeof(real)*(nr));
	weight = new real[nr];
	memset(weight, 0, sizeof(real)*(nr));

	delta = new real[nr];
	memset(delta, 0, sizeof(real)*(nr));

	Sij = new real*[nr + 1];
	for (i = 0; i < (nr + 1); i++) {
		Sij[i] = new real[nr + 1];
		memset(Sij[i], 0, sizeof(real)*(nr + 1));
	}

	Sijk = new real**[nr + 1];
	for (i = 0; i < (nr + 1); i++) {
		Sijk[i] = new real*[nr + 1];
		for (j = 0; j < (nr + 1); j++) {
			Sijk[i][j] = new real[nr];
			memset(Sijk[i][j], 0, sizeof(real)*(nr));
		}
	}////

	Pij = new real*[nr];
	for (i = 0; i < nr; i++) {
		Pij[i] = new real[nr];
		memset(Pij[i], 0, sizeof(real)*nr);
	}

	vol = new real[nr];
	memset(vol, 0, sizeof(real)*(nr));
	ci = new real[nr];
	memset(ci, 0, sizeof(real)*(nr));

	rho = new real*[nr + 1];
	for (i = 0; i < (nr + 1); i++) {
		rho[i] = new real[ngauss];
		memset(rho[i], 0, sizeof(real)*(ngauss));
	}
	gamma = 0.;

	// variable for matrix calculation
	/*b_y = new real[nr];
	memset(b_y, 0, sizeof(real)*(nr));
	
	b_x = new real[nr];
	memset(b_x, 0, sizeof(real)*(nr));

	Yi = new real[nr];
	memset(Yi, 0, sizeof(real)*(nr));

	Xiktmp = new real[nr];
	memset(Xiktmp, 0, sizeof(real)*(nr));

	Ya = new real[nr];
	memset(Ya, 0, sizeof(real)*(nr));

	xk = new real[nr];
	memset(xk, 0, sizeof(real)*(nr));

	sourceFlux = new real[nr];
	memset(sourceFlux, 0, sizeof(real)*(nr));

	J_Flux = new real[nr];
	memset(J_Flux, 0, sizeof(real)*(nr));

	Flux = new real[nr];
	memset(Flux, 0, sizeof(real)*(nr));

	A = new real*[nr];
	for (i = 0; i < nr; i++) {
		A[i] = new real[nr];
		memset(A[i], 0, sizeof(real)*nr);
	}

	Xik = new real*[nr];
	for (i = 0; i < nr; i++) {
		Xik[i] = new real[nr];
		memset(Xik[i], 0, sizeof(real)*nr);
	}

	Xika = new real*[nr];
	for (i = 0; i < nr; i++) {
		Xika[i] = new real[nr];
		memset(Xika[i], 0, sizeof(real)*nr);
	}*/

	///
	removal = new real[nr];
	memset(sigs, 0, sizeof(real)*(nr));

}

real CollisionProbability::calBickeley(real x_in)
{
	//	FUNCTION KI3(XX)
	//		C / HEBE, 02 Jan 1996 ---- - By Rudi Stamm'ler        *        *       
	//		C   ACCURACY 1.0E-5, UP TO XX = 9.0.

	double KI3;
	double X = fabs(x_in);
	int i;

	if (X > 0.99999)
	{
		if (X > 9.0) {
			KI3 = 0.0;
		}
		else
		{
			i = floor(2.5*(X - 0.99998)) + 1;

			switch (i)
			{
				//	% C                           ** RANGE 1.0 - 1.4 **
			case 1:
				KI3 = ((-0.05337485*X + 0.3203223)*X - 0.7538355)*X + 0.7247294;
				break;
				//% C                           ** RANGE 1.4 - 1.8 **
			case  2:
				KI3 = ((-0.03146833*X + 0.2295280)*X - 0.6279752)*X + 0.6663720;
				break;
				//% C                           ** RANGE 1.8 - 2.2 **
			case 3:
				KI3 = ((-0.01906198*X + 0.1631667)*X - 0.5094124)*X + 0.5956163;
				break;
				//% C                           ** RANGE 2.2 - 2.6 **
			case 4:
				KI3 = ((-0.01174752*X + 0.1152418)*X - 0.4046007)*X + 0.5191031;
				break;
				//% C                           ** RANGE 2.6 - 3.0 **
			case 5:
				KI3 = ((-0.007328415*X + 0.08097913)*X - 0.3159648)*X + 0.4425954;
				break;
				//% C                           ** RANGE 3.0 - 3.4 **
			case 6:
				KI3 = ((-0.004617254*X + 0.05669960)*X - 0.2434341)*X + 0.3703178;
				break;
				//% C                           ** RANGE 3.4 - 3.8 **
			case 7:
				KI3 = (0.007923547*X - 0.07158569)*X + 0.1684022;
				break;
				//% C                           ** RANGE 3.8 - 4.2 **
			case 8:
				KI3 = (0.005095111*X - 0.05016344)*X + 0.1278307;
				break;
				//% C                           ** RANGE 4.2 - 4.6 **
			case 9:
				KI3 = (0.003286040*X - 0.03501524)*X + 0.09611422;
				break;
				//% C                           ** RANGE 4.6 - 5.0 **
			case 10:
				KI3 = (0.002126242*X - 0.02437465)*X + 0.07170491;
				break;
				//% C                           ** RANGE 5.0 - 5.8 **
			case 11:
				KI3 = (0.001123687*X - 0.01425519)*X + 0.04616317;
				break;
				//% C                           ** RANGE 5.8 - 6.6 **
			case 12:
				KI3 = (0.001123687*X - 0.01425519)*X + 0.04616317;
				break;
				//% C                           ** RANGE 5.8 - 6.6 **
			case 13:
				KI3 = (4.762937E-4*X - 6.810124E-3)*X + 0.02475115;
				break;
				//% C                           ** RANGE 6.6 - 7.4 **
			case 14:
				KI3 = (4.762937E-4*X - 6.810124E-3)*X + 0.02475115;
				break;
				//% C                           ** RANGE 6.6 - 7.4 **
			case 15:
				KI3 = (2.031843E-4*X - 3.232035E-3)*X + .01302864;
				break;
				//% C                           ** RANGE 7.4 - 8.2 **
			case 16:
				KI3 = (2.031843E-4*X - 3.232035E-3)*X + .01302864;
				break;
				//% C                           ** RANGE 7.4 - 8.2 **
			case 17:
				KI3 = (8.701440E-5*X - 1.524126E-3)*X + 6.749972E-3;
				break;
				//% C                           ** RANGE 8.2 - 9.0 **
			case 18:
				KI3 = (8.701440E-5*X - 1.524126E-3)*X + 6.749972E-3;
				break;
				//% C                           ** RANGE 8.2 - 9.0 **
			case 19:
				KI3 = (3.742673E-5*X - 7.157367E-4)*X + 3.454768E-3;
				break;
			case 20:
				KI3 = (3.742673E-5*X - 7.157367E-4)*X + 3.454768E-3;
				break;
			case 21:
				KI3 = (3.742673E-5*X - 7.157367E-4)*X + 3.454768E-3;
				break;
			}
		}
	}
	else {
		i = floor(20.0*X + 1.00001);
		switch (i)
		{
			//	% C                          ** RANGE 0.00 - 0.05 **
		case 1:
			KI3 = (0.7266088*X - 0.9990226)*X + 0.7853961;
			break;
			//% C                          ** RANGE 0.05 - 0.10 **
		case 2:
			KI3 = (0.6466375*X - 0.9912340)*X + 0.7852024;
			break;
			//% C                          ** RANGE 0.10 - 0.15 **
		case 3:
			KI3 = (0.5856605*X - 0.9791293)*X + 0.7845986;
			break;
			//% C                          ** RANGE 0.15 - 0.20 **
		case 4:
			KI3 = (0.5346648*X - 0.9638914)*X + 0.7834577;
			break;
			//% C                          ** RANGE 0.20 - 0.25 **
		case 5:
			KI3 = (0.4907827*X - 0.9463843)*X + 0.7817094;
			break;
			//% C                          ** RANGE 0.25 - 0.30 **
		case 6:
			KI3 = (0.4521752*X - 0.9271152)*X + 0.7793031;
			break;
			//% C                          ** RANGE 0.30 - 0.35 **
		case 7:
			KI3 = (0.4177388*X - 0.9064822)*X + 0.7762107;
			break;
			//% C                          ** RANGE 0.35 - 0.40 **
		case 8:
			KI3 = (0.3869945*X - 0.8849865)*X + 0.7724519;
			break;
			//% C                          ** RANGE 0.40 - 0.45 **
		case 9:
			KI3 = (0.3590753*X - 0.8626685)*X + 0.7679903;
			break;
			//% C                          ** RANGE 0.45 - 0.50 **
		case 10:
			KI3 = (0.3338676*X - 0.8400133)*X + 0.7628988;
			break;
			//% C                          ** RANGE 0.50 - 0.60 **
		case 11:
			KI3 = (0.2998569*X - 0.8054172)*X + 0.7540982;
			break;
			//% C                          ** RANGE 0.60 - 0.70 **
		case 12:
			KI3 = (0.2998569*X - 0.8054172)*X + 0.7540982;
			break;
			//% C                          ** RANGE 0.60 - 0.70 **
		case 13:
			KI3 = (0.2609154*X - 0.7587821)*X + 0.7401279;
			break;
			//% C                          ** RANGE 0.70 - 0.80 **
		case 14:
			KI3 = (0.2609154*X - 0.7587821)*X + 0.7401279;
			break;
			//% C                          ** RANGE 0.70 - 0.80 **
		case 15:
			KI3 = (0.2278226*X - 0.7125290)*X + 0.7239594;
			break;
			//% C                          ** RANGE 0.80 - 0.90 **
		case 16:
			KI3 = (0.2278226*X - 0.7125290)*X + 0.7239594;
			break;
			//% C                          ** RANGE 0.80 - 0.90 **
		case 17:
			KI3 = (0.1994999*X - 0.6672761)*X + 0.7058777;
			break;
			//% C                          ** RANGE 0.90 - 1.00 **
		case 18:
			KI3 = (0.1994999*X - 0.6672761)*X + 0.7058777;
			break;
			//% C                          ** RANGE 0.90 - 1.00 **
		case 19:
			KI3 = (0.1751248*X - 0.6234536)*X + 0.6861762;
			break;
		case 20:
			KI3 = (0.1751248*X - 0.6234536)*X + 0.6861762;
			break;
		}
	}

	return KI3;
}

void CollisionProbability::gaussJacobiQuad(int inGaussNum, real * x_in, real *weight_in)
{	//Gauss Jacobi Quadrature

	switch (inGaussNum)
	{
	case 1:
		x_in[0] = 0.6666666667;

		weight_in[0] = 0.5000000000;

		break;
	case 2:
		x_in[0] = 0.3550510257;
		x_in[1] = 0.8449489743;

		weight_in[0] = 0.1819586183;
		weight_in[1] = 0.3180413817;

		break;
	case 3:
		x_in[0] = 0.2123405382;
		x_in[1] = 0.5905331356;
		x_in[2] = 0.9114120405;

		weight_in[0] = 0.0698269799;
		weight_in[1] = 0.2292411064;
		weight_in[2] = 0.2009319137;

		break;
	case 4:
		x_in[0] = 0.1397598643;
		x_in[1] = 0.4164095676;
		x_in[2] = 0.7231569864;
		x_in[3] = 0.9428958039;

		weight_in[0] = 0.0311809710;
		weight_in[1] = 0.1298475476;
		weight_in[2] = 0.2034645680;
		weight_in[3] = 0.1355069134;

		break;
	case 5:
		x_in[0] = 0.0985350858;
		x_in[1] = 0.3045357266;
		x_in[2] = 0.5620251898;
		x_in[3] = 0.8019865821;
		x_in[4] = 0.9601901429;

		weight_in[0] = 0.0157479145;
		weight_in[1] = 0.0739088701;
		weight_in[2] = 0.1463869871;
		weight_in[3] = 0.1671746381;
		weight_in[4] = 0.0967815902;

		break;
	default:
		cout << "check the ngauss in the input file is proper." << endl;
		break;
	}
}

void CollisionProbability::gaussianQuad(int inGaussNum, real * x_in, real * weight_in)
{
	//Gauss Quadrature

	switch (inGaussNum)
	{
	case 2:
		x_in[0] = 0.577350269189626;
		x_in[1] = -x_in[0];

		weight_in[0] = 1.0;
		weight_in[1] = weight_in[0];

		break;
	case 3:
		x_in[0] = 0.774596668241483;
		x_in[1] = 0.;
		x_in[2] = -x_in[0];

		weight_in[0] = 0.555555555555556;
		weight_in[1] = 0.888888888888889;
		weight_in[2] = weight_in[0];

		break;
	case 4:
		x_in[0] = 0.339981043584856;
		x_in[1] = 0.861136311594053;
		x_in[2] = -x_in[0];
		x_in[3] = -x_in[1];

		weight_in[0] = 0.652145154862546;
		weight_in[1] = 0.347854845137454;
		weight_in[2] = weight_in[0];
		weight_in[3] = weight_in[1];

		break;
	case 5:
		x_in[0] = 0.538469310105683;
		x_in[1] = 0.906179845938664;
		x_in[2] = 0.;
		x_in[3] = -x_in[0];
		x_in[4] = -x_in[1];

		weight_in[0] = 0.478628670499366;
		weight_in[1] = 0.236826885056189;
		weight_in[2] = 0.568888888888889;
		weight_in[3] = weight_in[0];
		weight_in[4] = weight_in[1];

		break;
	default:
		cout << "check the ngauss in the input file is proper." << endl;
		break;
	}

}

void CollisionProbability::calCPKernels()
{
	int i, j, k,m;
	real temp, sum;
	real y;
	
	delta[0] = rad[0];

	for (i = 0; i < nr; i++) {
		sigr[i] = sigt[i] - sigs[i ];
	}
	for (i = 1; i < nr; i++) {
		delta[i] = rad[i] - rad[i - 1];
	}
	vol[0] = rad[0] * rad[0] * _pi;
	for (i = 1; i < nr; i++) {
		vol[i] = (rad[i] * rad[i] - rad[i - 1] * rad[i - 1])*_pi;
	}
	for (i = 0; i < nr; i++) {
		ci[i] = (sigs[i] / sigt[i]);
	}

	for (k = 0; k < nr; k++) {
		// calculate rho
		for (m = 0; m < ngauss; m++) {
			for (i = k; i < nr; i++) {
				if (i == k) {
					y = rad[k] - x_i[m] * x_i[m] * delta[k];
					rho[i + 1][m] = (sqrt(rad[i] * rad[i] - y*y))*sigt[i];// +rho[i][m];
				}
				else {
					y = rad[k] - x_i[m] * x_i[m] * delta[k];
					rho[i + 1][m] = (sqrt(rad[i] * rad[i] - y*y) - sqrt(rad[i - 1] * rad[i - 1] - y*y))*sigt[i] + rho[i][m];
				}
			}
		}

		//cal Sijk
		for (i = k + 1; i < (nr + 1); i++) {
			for (j = i; j < (nr + 1); j++) {
				for (m = 0; m < ngauss; m++) {
					Sijk[i][j][k] += 2 * delta[k] * (weight[m])*(calBickeley(rho[j][m]+ rho[i][m]) - calBickeley(rho[j][m]- rho[i][m]));
				}
				Sij[i][j] += Sijk[i][j][k];
			}
		}
	}

	for (i = 1; i < (nr + 1); i++) {
		for (j = 0; j < i; j++) {
			Sij[i][j] = Sij[j][i];
		}
	}

	for (i = 0; i < nr; i++) {
		for (j = 0; j < nr; j++) {
			//calPij(i, j);
			if (i == j) sum = vol[i] * sigt[i];
			else
			{
				sum = 0.;
			}
			Pij[i][j] = sum + 2 * (Sij[i + 1][j + 1] + Sij[i][j] - Sij[i][j + 1] - Sij[i + 1][j]);
		}
	}
}

void CollisionProbability::calCPKernels_GaussianQuad()
{
	int i, j, k, m;
	real temp, sum;
	real y;

	delta[0] = rad[0];

	for (i = 0; i < nr; i++) {
		sigr[i] = sigt[i] - sigs[i];
	}
	for (i = 1; i < nr; i++) {
		delta[i] = rad[i] - rad[i - 1];
	}
	vol[0] = rad[0] * rad[0] * _pi;
	for (i = 1; i < nr; i++) {
		vol[i] = (rad[i] * rad[i] - rad[i - 1] * rad[i - 1])*_pi;
	}
	for (i = 0; i < nr; i++) {
		ci[i] = (sigs[i] / sigt[i]);
	}

	for (k = 0; k < nr; k++) {
		// calculate rho
		for (m = 0; m < ngauss; m++) {
			for (i = k; i < nr; i++) {
				if (i == k) {
					y = rad[k] - x_i[m] * x_i[m] * delta[k];
					//y = rad[k] - x_i[m] * x_i[m] * delta[k];
					rho[i + 1][m] = (sqrt(rad[i] * rad[i] - y*y))*sigt[i];// +rho[i][m];
				}
				else {
					y = rad[k] - x_i[m] * x_i[m] * delta[k];
					rho[i + 1][m] = (sqrt(rad[i] * rad[i] - y*y) - sqrt(rad[i - 1] * rad[i - 1] - y*y))*sigt[i] + rho[i][m];
				}
			}
		}

		//cal Sijk
		for (i = k + 1; i < (nr + 1); i++) {
			for (j = i; j < (nr + 1); j++) {
				for (m = 0; m < ngauss; m++) {
					Sijk[i][j][k] += 1./2. * delta[k] * (weight[m])*(calBickeley(rho[j][m] + rho[i][m]) - calBickeley(rho[j][m] - rho[i][m]));
				}
				Sij[i][j] += Sijk[i][j][k];
			}
		}
	}

	for (i = 1; i < (nr + 1); i++) {
		for (j = 0; j < i; j++) {
			Sij[i][j] = Sij[j][i];
		}
	}

	for (i = 0; i < nr; i++) {
		for (j = 0; j < nr; j++) {
			//calPij(i, j);
			if (i == j) sum = vol[i] * sigt[i];
			else
			{
				sum = 0.;
			}
			Pij[i][j] = sum + 2 * (Sij[i + 1][j + 1] + Sij[i][j] - Sij[i][j + 1] - Sij[i + 1][j]);
		}
	}
}

void CollisionProbability::solveSystem()
{

	int i, j, k;

	b_y = new real[nr];
	memset(b_y, 0, sizeof(real)*(nr));

	b_x = new real[nr];
	memset(b_x, 0, sizeof(real)*(nr));

	Yi = new real[nr];
	memset(Yi, 0, sizeof(real)*(nr));

	Xiktmp = new real[nr];
	memset(Xiktmp, 0, sizeof(real)*(nr));

	Ya = new real[nr];
	memset(Ya, 0, sizeof(real)*(nr));

	xk = new real[nr];
	memset(xk, 0, sizeof(real)*(nr));

	sourceFlux = new real[nr];
	memset(sourceFlux, 0, sizeof(real)*(nr));

	J_Flux = new real[nr];
	memset(J_Flux, 0, sizeof(real)*(nr));

	Flux = new real[nr];
	memset(Flux, 0, sizeof(real)*(nr));

	A = new real*[nr];
	for (i = 0; i < nr; i++) {
	A[i] = new real[nr];
	memset(A[i], 0, sizeof(real)*nr);
	}

	Xik = new real*[nr];
	for (i = 0; i < nr; i++) {
	Xik[i] = new real[nr];
	memset(Xik[i], 0, sizeof(real)*nr);
	}

	Xika = new real*[nr];
	for (i = 0; i < nr; i++) {
	Xika[i] = new real[nr];
	memset(Xika[i], 0, sizeof(real)*nr);
	}


	real totalx;
	real Jout, Jin;

	real deno;
	real temp;
	
	// construct matrix A

	for (i = 0; i < nr; i++) {
		for (j = 0; j < nr; j++) {
			A[i][j] = -Pij[i][j] * ci[j];
		}
		A[i][i] = A[i][i] + sigt[i] * vol[i];
	}

	// construct vector b_y for Yi
	for (i = 0; i < nr; i++) {
		temp = sigt[i] * vol[i];
		for (j = 0; j < nr; j++) {
			temp = temp -Pij[i][j];
		}
		b_y[i] = 4.0 / (2*_pi*rad[nr-1])*temp;
	}

	// LU Factorization
	luFactorization(A);

	//calculate Yi
	solveLU(A, Yi, b_y);


	// construct vector b_x for Xik

	for (k = 0; k < nr; k++) {
		for (i = 0; i < nr; i++) {
			b_x[i] = Pij[k][i] / sigt[k];
		}
		//calculate Xiktmp
		solveLU(A, Xiktmp, b_x);

		for (i = 0; i < nr; i++) {
			Xik[i][k] = Xiktmp[i];
		}
	}
	
	//cal gamma
	for (i = 0; i < nr; i++) {
		gamma += (sigt[i] - sigs[i])*Yi[i]*vol[i];
	}

	deno = 1 - (1 - gamma)*albedo;
	for (i = 0; i < nr; i++) {
		Ya[i] = Yi[i] / deno;
	}

	for (i = 0; i < nr; i++) {
		xk[i] = _pi*rad[nr - 1] / 2 * vol[i] * Yi[i] ;
	}

	totalx = 0;
	for (i = 0; i < nr; i++) {
		totalx += qi[i] * xk[i] ;
	}

	for (k = 0; k < nr; k++) {
	
		for (i = 0; i < nr; i++) {
			Xika[i][k] = Xik[i][k]+albedo*xk[k] *Ya[i] ;
		}
	}

	Jout = (totalx + (1 - gamma)*jext) / (1 - albedo*(1 - gamma));
	Jin = (albedo*totalx + jext) / (1 - albedo*(1 - gamma));

	for (i = 0; i < nr; i++) {

		for (k = 0; k < nr; k++) {
			sourceFlux[i]  += Xika[i][k] * qi[k];
		}
		J_Flux[i]  = Ya[i]  * jext ;
		Flux[i]  = J_Flux[i]  + sourceFlux[i] ;
	}
	

	//////////////////////////////////////////////////////////////////////////////////////////////
/*
	ofstream fout("variable4debug.txt");
	


	//fout.setf(ios::fixed);
	fout.setf(ios::scientific);
	fout << setprecision(4);

	//파일에 결과값을 입력한다.
	fout << "Sijk" << endl;
	for (k = 0; k < nr; k++) {
		fout << "k : " << k << endl;
		for (i = 0; i < nr + 1; i++) {
			for (j = 0; j < nr + 1; j++) {
				fout << Sijk[i][j][k] << "\t";
			}
			fout << endl;
		}
		fout << endl;
	}
	fout << endl;

	//파일에 결과값을 입력한다.
	fout << "Sij" << endl;
	for (i = 0; i < nr + 1; i++) {
		for (j = 0; j < nr + 1; j++) {
			fout << Sij[i][j] << "\t";
		}
		fout << endl;
	}
	fout << endl;

	fout << "Pij" << endl;
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nr; j++) {
			fout << Pij[i][j] << "\t";
		}
		fout << endl;
	}
	fout << endl;

	fout << "Yi" << endl;
	for (i = 0; i < nr; i++) {
		fout << Yi[i]  << "\t";
	}
	fout << endl<<endl;

	fout << "Ya" << endl;
	for (i = 0; i < nr; i++) {
		fout << Ya[i]  << "\t";
	}
	fout << endl << endl;

	fout << "A" << endl;
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nr; j++) {
			fout << A[i][j]<< "\t";
			
		}
		fout << endl;
	}
	fout << endl;

	fout << "Xik" << endl;
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nr; j++) {
			fout << Xik[i][j] << "\t";

		}
		fout << endl;
	}
	fout << endl;

	fout << "Xika" << endl;
	for (i = 0; i < nr; i++) {
		for (j = 0; j < nr; j++) {
			fout << Xika[i][j] << "\t";

		}
		fout << endl;
	}
	fout << endl;

	fout << "xk" << endl;
	for (i = 0; i < nr; i++) {
		fout << xk[i]  << "\t";
	}
	fout << endl << endl;

	fout << "totalx" << endl;
	fout << totalx<< "\t";
	fout << endl;

	fout << "Jout" << endl;
	fout << Jout << "\t";
	fout << endl;

	fout << "Jin" << endl;
	fout << Jin << "\t";
	fout << endl;

	fout << "gamma" << endl;
	fout << gamma << "\t";
	fout << endl;
	
	fout << "Flux" << endl;
	for (i = 0; i < nr; i++) {
		fout << sourceFlux[i]  << "\t" <<J_Flux[i]  << "\t" <<Flux[i]  << endl;
	}
	fout << endl << endl;
*/

	//Output.txt <- Result
	
	real Rtot=0;
	real Jnet;
	Jnet = Jout - Jin;

	real qtot = 0;

	for (i = 0; i < nr; i++) {
		removal[i]  = sigr[i] * vol[i] * Flux[i] ;
		Rtot += removal[i] ;
		qtot += vol[i] * qi[i];
	}
	
	string outName;
	outName = "Output_of_";
	outName.append(inName);
	outName.append(".txt");

	ofstream fout1(outName);

	int width,precision;
	precision = 6;
	width = precision+10;
	//fout1.setf(ios::fixed);
	fout1.setf(ios::scientific);
	fout1 << setprecision(precision);
	fout1 << "Albedo =  " << albedo << "     Jext =   " << jext << endl;
	fout1 << "     Jnet =" <<setw(width)<< Jnet<< "     jout =" <<setw(width)<< Jout << "     jin =" << setw(width) << Jin << endl;
	fout1 << "      Region      Radius      SigRem     SigScat      Source    SRC Flux      J Flux  Total Flux     Removal" << endl;
	for (i = 0; i < nr; i++)
		fout1 << setw(width)<<(i+1)<<setw(width)<< rad[i] << setw(width) << sigr[i] << setw(width) << sigs[i] << setw(width) << qi[i]<< setw(width)<<sourceFlux[i] << setw(width) << J_Flux[i]  << setw(width) << Flux[i] <<setw(width) <<removal[i]  << endl;
	fout1 << "     Rtot =" << setw(14) << Rtot << "     jnet =" << setw(14) << Jnet << " sum =" << setw(14) << qtot << " qtot "<< qtot <<endl;

	fout1.close();
}


void CollisionProbability::luFactorization(real ** matA)
{
	int i, j, k, i_temp, j_temp;
	real temp;

	for (j = 0; j < nr-1 ; j++) {		// 0:nr-2 -> 1:nr-1
		for (i = j+1; i < nr; i++) {		// j:nr-1 -> j+1:nr
			matA[i][j] = matA[i][j] / matA[j][j];
		}
		for (k = j+1; k < nr; k++) {		// j:nr-1 -> j+1:nr
			for (i = j+1; i < nr; i++) {	// j:nr-1 -> j+1:nr
				matA[k][i] = matA[k][i] - matA[k][j] * matA[j][i];
			}
		}
	}

	//for j=1:n-1
	//    for i=j+1:n
	//        A(i,j)=A(i,j)/A(j,j);
	//    end
	//    for k=j+1:n
	//        for i=j+1:n
	//            A(k,i)=A(k,i)-A(k,j)*A(j,i);
	//        end
	//    end
	//end
}

void CollisionProbability::solveLU(real ** matA, real * vec_x, real * vec_b)
{
	int i, j, k, i_temp,j_temp;
	real temp;

	// Forward substitution
	vec_x[0] = vec_b[0];
	for (i = 1; i < nr;i++) {
		temp = 0;
		for (j = 0; j < i ;j++) {
			temp = temp + matA[i][j]*vec_x[j];
		}
		vec_x[i] = (vec_b[i] - temp);
	}
	// Backward substitution
	vec_x[nr-1] = vec_x[nr-1] / matA[nr-1][nr-1];
	for (i = nr-2; i>=0 ; i--) {
			temp = 0;
			for (j= nr-1; j >i ; j--) {
					temp = temp + matA[i][j]*vec_x[j];
			}
		vec_x[i] = (vec_x[i] - temp) / matA[i][i];
	}
}


