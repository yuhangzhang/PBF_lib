
#include "engine.h"
#include "graph.h"
#include "QPBpoly.h"

typedef struct cQPBFlist
{
	QPBpoly* A_i;
	int cid;
	cQPBFlist* next;
}cQPBFlist;


class Learning_QPBF
{
public:
	Learning_QPBF();
	int add_cQPBF(QPBpoly* A_i); //return the id of the added component
	void remove_cQPBF(int cid);//remove component-cid
	void learn(Matrix<bool,Dynamic,1> y, Matrix<double,Dynamic,1> coeff, Matrix<double,Dynamic,1> lambda);
	int numvar();
	int numcomp();
	double optimize(Matrix<bool,Dynamic,1>& y);
	void para(Matrix<double,Dynamic,1>& p);
private:
	Matrix<double,Dynamic,1> _para;
	cQPBFlist* _componentlist;
	int _numvar;
};



