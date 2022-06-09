#pragma once
#include "define.h"

class Matrix;
class StripeMatrix;

// Iterative methods to solve Ax = b
namespace IterativeLinearSolver
{
	const double ConvCri = 1e-8;
	void Jacobi(const StripeMatrix& A, const vector<double>& b, vector<double>& x);
	void GaussSeidal(const StripeMatrix& A, const vector<double>& b, vector<double>& x);
	void SOR(const StripeMatrix& A, const vector<double>& b, double w, vector<double>& x);
}

static void PrintError(const string func, const string msg = "", ostream& out = cout) {
	out << "##### Error occurs at  [" << func << "]" << endl;
	if (!msg.empty())
		out << "##### [Error Message] : " << msg << endl;
}

class Matrix
{
protected:
	//----------------------------------------------------------------------------------------------------------------
	// Member Variables
	//----------------------------------------------------------------------------------------------------------------
	vector<vector<double>> val;
	unsigned nRow, nCol;

public:
	//----------------------------------------------------------------------------------------------------------------
	// Constructor & Destructor
	//----------------------------------------------------------------------------------------------------------------			
	Matrix();
	Matrix(int n);
	Matrix(int n, int m);
	Matrix(const Matrix& origin);
	virtual ~Matrix();
	//----------------------------------------------------------------------------------------------------------------
	// Operators
	//----------------------------------------------------------------------------------------------------------------
	const double& operator()(int row, int col) const { return val[row][col]; }
	double& operator()(int row, int col) { return const_cast<double&>(static_cast<const Matrix&>(*this)(row, col)); }
	const vector<double>& operator[](int row) const { return val[row]; }
	vector<double>& operator[](int row) { return const_cast<vector<double>&>(static_cast<const Matrix&>(*this)[row]); }
	const Matrix& operator = (const vector<double>& rowVec);
	const Matrix& operator = (const Matrix& origin);
	Matrix operator + (const Matrix& origin) const;
	Matrix operator - (const Matrix& origin) const;
	Matrix operator * (const Matrix& origin) const;
	Matrix operator * (double s) const;
	friend Matrix operator * (double s, const Matrix& origin);
	//----------------------------------------------------------------------------------------------------------------
	// Member Access Functions
	//----------------------------------------------------------------------------------------------------------------
	double GetVal(int i, int j) const { return val[i][j]; }
	unsigned GetNRow() const { return nRow; }
	unsigned GetNCol() const { return nCol; }
	//----------------------------------------------------------------------------------------------------------------
	// 
	//----------------------------------------------------------------------------------------------------------------
	void AllocateMemory(int nRow, int nCol);
	friend static void CheckDimension(int n1, int m1, int n2, int m2, char op = '=');
	friend static void CheckDimension(const Matrix& A, const Matrix& B, char op = '=');
	void Transpose();
	void PrintMatrix(ostream& out, int precision = 3, char type = 'f') const;
	//-----------------------------------------------------------------------------------------------------------------
	//Methods
	//-----------------------------------------------------------------------------------------------------------------
	void LU(Matrix& L, Matrix& U, Matrix& P, int& d) const;
	void LU(Matrix& L, Matrix& U, Matrix& P) const { int d; LU(L, U, P, d); }
	void LU(Matrix& L, Matrix& U) const { Matrix dum(nRow,nCol); int d; LU(L, U, dum, d); }
	double Det() const;
	void LinearSolver(Matrix&X, Matrix&b) const; //// A.LinearSolver(X,b);
	void Inv(Matrix& X) const;

private:
	bool OutOfRange(int i, int j) const { return (i < 0 || (int)nRow <= i || j < 0 || (int)nCol <= j) ? true : false; }

};

class StripeMatrix :public Matrix
{
private:
	vector<vector<int>> interval; // interval[i][j] : interval from diagonal eleme ith row jth column element

public:
	StripeMatrix(int nRow, int nCol);
	StripeMatrix(int N, const vector<int>& interval_);
	virtual ~StripeMatrix() {};
	void SetRowInter(int iRow, const vector<int>& rowInter) { interval[iRow] = rowInter; }
	void SetRowVec(int iRow, const vector<double>& rowVec) { val[iRow] = rowVec; }
	void SetColVec(int iCol, const vector<double>& colVec);
	virtual const StripeMatrix& operator=(const StripeMatrix& origin);
	void Multiply(const vector<double>& x, vector<double>& b) const;
	virtual void Transpose();
	friend void IterativeLinearSolver::Jacobi(const StripeMatrix& A, const vector<double>& b, vector<double>& x);
	friend void IterativeLinearSolver::GaussSeidal(const StripeMatrix& A, const vector<double>& b, vector<double>& x);
	friend void IterativeLinearSolver::SOR(const StripeMatrix& A, const vector<double>& b, double w, vector<double>& x);

private:
};