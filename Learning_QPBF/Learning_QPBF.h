#include "QPBpoly.h"
#include "engine.h"


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
	void add_cQPBF(QPBpoly A_i);
	void learn(Matrix<bool,Dynamic,1> y, Matrix<double,Dynamic,1> w);
private:
	Matrix<double,Dynamic,1> _para;
	cQPBFlist* _componentlist;
	int _numvar;
};



