#include <map>
#include <Eigen/core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
using namespace Eigen;
/*
Every object of QPBpoly stores a
Quadratic pseudo-Boolean function in a polynomial form. 

After initialization, the QPBF is an empty function. 

Modify the linear and bilinear terms in a QPBF with addTerm1 and addTerm2.


*/

typedef std::pair<std::pair<int,int>,double> QPBTerm;
typedef std::map<std::pair<int,int>,double> QPBF;

class QPBpoly
{
public:
  	QPBpoly(int numvar);//numvar is the number of boolean variables
	void addTerm1(int v0, double coeff);
	void addTerm2(int v0, int v1, double coeff);
	double getTerm1(int v0);
	double getTerm2(int v0, int v1);
	
	bool isSubmodular();
	void delTerm1(int v0);//set coefficient of term x_i to 0
	void delTerm2(int v0, int v1);

	QPBF::iterator firstTerm();
	QPBF::iterator lastTerm();//not really a term, but a node behind the last

	QPBpoly operator+(QPBpoly qpbf);
	QPBpoly operator-(QPBpoly qpbf);
	QPBpoly operator*(double u);
	double& operator()(int v0, int v1);//use () operators only if the matrix is not sparse or large
	double& operator()(int v0);

	void clear();//remove all terms
	void clean();//remove 0 terms

	int size();//number of nonzero terms

	int numvar();//return the number of boolean variables
	int updatenumvar(int numvar);

	double evaluate(Matrix<bool,Dynamic,1> y);

private:
	int _numvar;
	//_QPBcoeff stores the coefficient matrix A of the polynomial form
	//e.g. f(x)=x'Ax
	//It is a map from the row and column to the coefficient
	//we use a pair<int,int> to couple the row and column
	QPBF _QPBcoeff;
};


