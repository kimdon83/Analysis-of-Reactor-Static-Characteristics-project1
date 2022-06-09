#include "Matrix.h"

Matrix::Matrix() :nCol(0), nRow(0)
{
}

Matrix::Matrix(int n)
{
	if (n < 0)
	{
		PrintError(__FUNCSIG__, "n < 0");
		exit(-1);
	}
	AllocateMemory(n, n);
}

Matrix::Matrix(int n, int m)
{
	if (n < 0 || m < 0)
	{
		PrintError(__FUNCSIG__, " n < 0 or m < 0");
		exit(-1);
	}
	AllocateMemory(n, m);
}

Matrix::Matrix(const Matrix& origin)
{
	nRow = origin.nRow; nCol = origin.nCol;
	val = origin.val;
}

Matrix::~Matrix()
{
}

const Matrix& Matrix::operator = (const vector<double>& rowVec)
{
	int m = rowVec.size();
	// dimension check
	if (!val.empty())
		CheckDimension(nRow, nCol, 1, m);
	else
		val.resize(1);

	nRow = 1; nCol = m;
	val[0] = rowVec;
	return *this;
}

const Matrix& Matrix::operator = (const Matrix& origin)
{
	// dimension check
	if (!val.empty())
		CheckDimension(*this, origin);

	nRow = origin.nRow; nCol = origin.nCol;
	val = origin.val;
	return *this;
}

Matrix Matrix::operator + (const Matrix& origin) const
{
	// dimension check
	CheckDimension(*this, origin);

	Matrix tmp(nRow, nCol);
	for (unsigned row = 0; row < nRow; row++)
	{
		for (unsigned col = 0; col < nCol; col++)
		{
			tmp[row][col] = (*this)[row][col] + origin[row][col];
		}
	}
	return tmp;
}

Matrix Matrix::operator - (const Matrix& origin) const
{
	// dimension check
	CheckDimension(*this, origin);

	Matrix tmp(nRow, nCol);
	for (unsigned row = 0; row < nRow; row++)
	{
		for (unsigned col = 0; col < nCol; col++)
		{
			tmp[row][col] = (*this)[row][col] - origin[row][col];
		}
	}
	return tmp;
}

Matrix Matrix::operator * (const Matrix& origin) const
{
	// dimension check
	CheckDimension(*this, origin, '*');

	Matrix tmp(nRow, origin.nCol);
	for (unsigned row = 0; row < nRow; row++)
	{
		for (unsigned col = 0; col < origin.nCol; col++)
		{

			for (unsigned k = 0; k < nCol; k++)
			{
				tmp[row][col] += (*this)[row][k] * origin[k][col];
			}
		}
	}
	return tmp;
}

Matrix Matrix:: operator * (double s) const
{
	Matrix tmp(nRow, nCol);
	for (unsigned row = 0; row < nRow; row++)
	{
		for (unsigned col = 0; col < nCol; col++)
		{
			tmp.val[row][col] = (*this)[row][col] * s;
		}
	}
	return tmp;
}

Matrix operator * (double s, const Matrix& origin)
{
	return origin * s;
}

void Matrix::AllocateMemory(int arg_nRow, int arg_nCol)
{
	nRow = arg_nRow; nCol = arg_nCol;
	val.resize(nRow);
	for (size_t row = 0; row < nRow; row++)
	{
		val[row].resize(nCol);
	}
}

static void CheckDimension(int n1, int m1, int n2, int m2, char op)
{
	switch (op)
	{
	case '=':
	case '+':
	case '-':
		if ((n1 == n2) && (m1 == m2))
			return;
		break;
	case '*':
		if ((m1 == n2))
			return;
		break;
	default:
		PrintError(__FUNCSIG__, "Current operator is not '=''+''-''*'");
		break;
	}
	string errMsg = "Dimension mismatch for operator ";
	errMsg.push_back(op);
	PrintError(__FUNCSIG__, errMsg);
	exit(-1);
	return;
}

void CheckDimension(const Matrix& A, const Matrix& B, char op)
{
	CheckDimension(A.nRow, A.nCol, B.nRow, B.nCol, op);
	return;
}

void Matrix::Transpose()
{
	Matrix tmp = (*this);
	if (nRow != nCol)
	{
		val.clear();
		swap(nRow, nCol);
		AllocateMemory(nRow, nCol);
	}

	for (unsigned row = 0; row < nRow; row++)
	{
		for (unsigned col = 0; col < nCol; col++)
		{
			(*this)[row][col] = tmp[col][row];
		}
	}
}

void Matrix::PrintMatrix(ostream& out, int precision, char type) const
{
	out.precision(precision);
	int width = precision + 4;
	switch (type)
	{
	case 'f':
	case 'F':
		out << fixed;
		break;
	case 'e':
	case 'E':
		out << scientific;
		out.setf(ios_base::uppercase);
		break;
	default:
		PrintError(__FUNCSIG__, "f/F : float , e/E : scientific");
		break;
	}
	for (unsigned row = 0; row < nRow; row++)
	{
		for (unsigned col = 0; col < nCol; col++)
		{
			out << setw(width) << (*this)[row][col];
		}
		out << endl;
	}
}

void Matrix::LU(Matrix& L, Matrix& U, Matrix& P, int& d) const
{
	//demension check
	if (!((nCol) == (nRow)) || nCol == 0)
	{
		PrintError(__FUNCSIG__, "not square Matrix!");
		exit(-1);
	}
	//initializing
	Matrix tmpU = *this;
	Matrix tmpL(nRow, nCol), tmpP(nRow, nCol), pivotP(nRow, nCol), eye(nRow, nCol);
	Matrix tmp(1, nRow);
	double mp, mi;
	double multiplier;
	unsigned ip;
	for (unsigned i = 0; i<nRow; i++)
	{
		tmpL[i][i] = 1;
		tmpP[i][i] = 1;
		pivotP[i][i] = 1;
		eye[i][i] = 1;
	}
	//main loop
	//pivoting
	int counter = 1; //pivot flag
	for (unsigned i = 0; i<nRow; i++)
	{
		mp = fabs(tmpU[i][i]);
		ip = i;
		for (unsigned j = i; j<nRow; j++)
		{
			mi = fabs(tmpU[j][i]);
			if (mp<mi)
			{
				ip = j;
				mp = mi;
			}
		}
		pivotP = eye;
		//need pivot?
		if (!(i == ip))
		{
			for (unsigned k = 0; k < nRow; k++)
			{
				tmp[0][k] = pivotP[i][k];
				pivotP[i][k] = pivotP[ip][k];
				pivotP[ip][k] = tmp[0][k];
			}
			counter *= (-1);
		}
		//P,A update
		tmpP = pivotP*tmpP;
		tmpU = pivotP*tmpU;
		//get L
		tmpL = (pivotP*tmpL)*pivotP;
		for (unsigned j = i + 1; j<nRow; j++)
		{
			tmpL[j][i] = ((tmpU[j][i]) / (tmpU[i][i]));
		}
		//get U
		for (unsigned j = i + 1; j<nRow; j++)
		{
			multiplier = ((tmpU[j][i]) / (tmpU[i][i]));
			for (unsigned k = 0; k<nRow; k++)
			{
				tmpU[j][k] = (tmpU[j][k]) - (multiplier*(tmpU[i][k]));
			}
		}
	}
	U = tmpU;
	L = tmpL;
	P = tmpP;
	d = counter;

}

void Matrix::LinearSolver(Matrix&X, Matrix&b) const
{
	if (!(nRow == b.nRow))
	{
		PrintError(__FUNCSIG__, "Dimension error!");
		exit(-1);
	}
	double sumod;
	Matrix x(nRow, 1), y(nRow, 1);
	Matrix tmpL(nRow, nCol), tmpU(nRow, nCol), tmpP(nRow, nCol);
	LU(tmpL, tmpU, tmpP);
	b = (tmpP*b);
	//Ly=b
	y[0][0] = b[0][0];
	for (unsigned i = 1; i < nRow; i++)
	{
		sumod = 0;
		for (unsigned j = 0; j < i; j++)
		{
			sumod += (tmpL[i][j] * y[j][0]);
		}
		y[i][0] = (b[i][0] - sumod) / tmpL[i][i];
	}
	//Ux=y
	x[nRow - 1][0] = y[nRow - 1][0] / tmpU[nRow - 1][nRow - 1];
	for (int i = (nRow - 2); i >= 0; i--)
	{
		sumod = 0;
		for (int j = nRow - 1; j >= i + 1; j--)
		{
			sumod += (tmpU[i][j] * x[j][0]);
		}
		x[i][0] = (y[i][0] - sumod) / tmpU[i][i];
	}

	X = x;
}

StripeMatrix::StripeMatrix(int nRow, int nCol)
{
	AllocateMemory(nRow, nCol);
	interval.resize(nRow);
	for (size_t i = 0; i < interval.size(); i++)
	{
		interval[i].resize(nCol);
	}
}

StripeMatrix::StripeMatrix(int N, const vector<int>& interval_)
:interval(vector<vector<int>>(N,interval_))
{
	AllocateMemory(N, interval_.size());
}

void StripeMatrix::SetColVec(int iCol, const vector<double>& colVec)
{
	for (unsigned row = 0; row < nRow; row++)
	{
		(*this)[row][iCol] = colVec[row];
	}
}

const StripeMatrix& StripeMatrix::operator = (const StripeMatrix& origin)
{
	// dimension check
	if (!val.empty())
		CheckDimension(*this, origin);

	nRow = origin.nRow; nCol = origin.nCol;
	val = origin.val; interval = origin.interval;
	return *this;
}

void StripeMatrix::Multiply(const vector<double>& x, vector<double>& b) const
{
	if (!(this->GetNRow() == x.size()))
	{
		PrintError(__FUNCSIG__, "Dimension error!");
		exit(-1);
	}
	b.resize(x.size());
	int fj; //j in full matrix

	for (unsigned i = 0; i < this->GetNRow(); i++)
	{
		b[i] = 0;
		for (unsigned j = 0; j < this->GetNCol(); j++)
		{
			fj = i + (*this).interval[i][j];
			b[i] += (*this)[i][j] * x[fj];
		}
	}
}

void StripeMatrix::Transpose() 
{
	StripeMatrix tmp = (*this);
	int fj;
	int nRow = this->nRow;
	int nCol = this->nCol;

	for (int j = 0; j < nCol; j++)
	{
		for (int i = 0; i < nRow; i++)
		{
			int interval = tmp.interval[i][j];
				fj = i + interval;
				//(*this)[i][j] = 0;
				//this->interval[i][j] = -9999999;
				if (interval == 0)
				{
					(*this)[i][j] = tmp[i][j];
				}
			if (!(interval==-9999999)&&!(interval==0))
			{
				
				if (nCol > 2) //L matrix
				{
					if (fj < nRow)
					{
						if (nCol == 8)//3D,2G
						{
							(*this)[fj][7 - j] = tmp[i][j];
							(*this).interval[fj][7 - j] = -interval;
						}
						if (nCol == 7)//3D,1G
						{
							(*this)[fj][6 - j] = tmp[i][j];
							(*this).interval[fj][6 - j] = -interval;
						}
					}
				}
				else //F matrix
				{
					(*this)[fj][1 - j] = tmp[i][j];
					(*this).interval[fj][1 - j] = -interval;
				}
				
			}
		}
	}
}

void IterativeLinearSolver::Jacobi(const StripeMatrix & A, const vector<double>& b, vector<double>& x)
{
	vector<double> Pre_x(x);
	double maxErr;
	double term, diag;
	int fj; // j in full matrix
	int n = 0;
	double delta = 0;
	double deltaSq = 0;
	int xSize = (int)x.size();
	double convCriSq = ConvCri*ConvCri;

	do
	{
		n += 1;
		Pre_x = x;
		deltaSq = 0;


		for (unsigned i = 0; i < A.GetNRow(); i++)
		{
			term = b[i]; diag = 0.;

			for (unsigned j = 0; j<A.GetNCol(); j++) {
				if (A.interval[i][j] == 0) diag = A[i][j];
				else {
					fj = i + A.interval[i][j];
					if (0 <= fj && fj < (int)x.size())
						term -= A[i][j] * Pre_x[fj];
				}

			}
			if (diag == 0.) {
				cout << "A diagonal term of the matrix is zero at CVector::SolveByPointJacobiIterativeMethod" << endl;
				exit(0);
			}
			x[i] = term / diag;
			delta = x[i] - Pre_x[i];
			deltaSq += delta*delta;
		}
		/*maxErr = 0;
		for (unsigned i = 0; i < A.GetNRow(); i++)
		{
		relErr = fabs(x[i] - Pre_x[i]) / Pre_x[i];
		if (relErr>maxErr) maxErr = relErr;
		}*/

		maxErr = deltaSq / inner_product(x.begin(), x.end(), x.begin(), 0.0);
	} while (maxErr > convCriSq);
	cout << " Inner Iter : " << n;
}

void IterativeLinearSolver::GaussSeidal(const StripeMatrix& A, const vector<double>& b, vector<double>& x)
{
	vector<double> Pre_x(x);
	double maxErr;
	double term, diag;
	int fj; // j in full matrix
	int n = 0;//conuter
	double delta = 0;
	double deltaSq = 0;
	int xSize = (int)x.size();
	double convCriSq = ConvCri*ConvCri;

	do
	{
		n += 1;
		Pre_x = x;
		deltaSq = 0;

		for (unsigned i = 0; i < A.GetNRow(); i++)
		{
			term = b[i];

			for (unsigned j = 0; j<A.GetNCol(); j++) {
				if (A.interval[i][j] == 0) diag = A[i][j];//diag = A[i][dim]
				else {
					fj = i + A.interval[i][j];
					if (0 <= fj && fj < xSize)
						term -= A[i][j] * x[fj];
				}	
			}
			x[i] = term / diag;
			delta = x[i] - Pre_x[i];
			deltaSq += delta*delta;
		}
		/*maxErr = 0;
		for (unsigned i = 0; i < A.GetNRow(); i++)
		{
			relErr = fabs(x[i] - Pre_x[i]) / Pre_x[i];
			if (relErr>maxErr) maxErr = relErr;
		}*/
		maxErr =deltaSq / inner_product(x.begin(), x.end(), x.begin(), 0.0);
	} while (maxErr > convCriSq);
	cout << " Inner Iter : " << n;
}

void IterativeLinearSolver::SOR(const StripeMatrix& A, const vector<double>& b, double w, vector<double>& x)
{
	vector<double> Pre_x(x);
	double maxErr;
	double term, diag;
	int fj; // j in full matrix
	int n = 0;//conuter
	double delta = 0;
	double deltaSq = 0;
	int xSize = (int)x.size();
	int nRow = A.GetNRow();
	int nCol = A.GetNCol();
	double convCriSq = ConvCri*ConvCri;

	do{
		n +=1;
		Pre_x = x;
		deltaSq = 0;


		for (int i = 0; i < nRow; i++)
		{
			term = b[i]; diag = 0.;
			for (int j = 0; j < nCol; j++)
			{
				if (A.interval[i][j] == 0) diag = A[i][j];//diag = A[i][dim]
				else {
					fj = i + A.interval[i][j];
					if (0 <= fj && fj < xSize)
						term -= A[i][j] * x[fj];
				}	
				x[i] = (1.0 - w)*Pre_x[i] + term*(w / diag);

			}
			delta = x[i] - Pre_x[i];
			deltaSq += delta*delta;
		}
		/*maxErr = 0;
		for (unsigned i = 0; i < A.GetNRow(); i++)
		{
		relErr = fabs(x[i] - Pre_x[i]) / Pre_x[i];
		if (relErr>maxErr) maxErr = relErr;
		}*/
		maxErr = deltaSq / inner_product(x.begin(), x.end(), x.begin(), 0.0);

	} while (maxErr > convCriSq);
	cout << " Inner Iter : " << n;
}