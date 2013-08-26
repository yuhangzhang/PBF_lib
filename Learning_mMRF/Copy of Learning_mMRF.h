

#include "./Learning_QPBF/Learning_QPBF.h"



class Learning_mMRF
{
public:
	Learning_mMRF(int nlabel);//initialized the number of labels
	int add_cQPBF(int label, QPBpoly* A_i); //return the id of the added component,
											//label specifies which label the added component function affects 
											//label<0 indicates all labels. 
	void remove_cQPBF(int cid);//remove component-cid, it only remove its effect, but do not change the numcomp
	void learn(Matrix<int,Dynamic,1> y);
	int numvar();
	int numcomp();
	double optimize(Matrix<int,Dynamic,1>& y);
private:
	Matrix<double,Dynamic,1> _para;
	QPBpoly _lctable;//tables recording the relations between labels and component functions
	cQPBFlist* _componentlist;
	int _numvar;
	int _numlabel;
};