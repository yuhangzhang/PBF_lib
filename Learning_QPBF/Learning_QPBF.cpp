#include "Learning_QPBF.h"


Learning_QPBF::Learning_QPBF()
{
	_componentlist = NULL;

	_para.resize(0,1);

	return;
}


void Learning_QPBF::add_cQPBF(QPBpoly A_i)
{
	cQPBFlist temp;
	temp.A_i = &A_i;
	temp.next= _componentlist;
	_componentlist = &temp;

	_para.resize(_para.rows()+1,1);

	return;
}



void Learning_QPBF::learn(Matrix<bool,Dynamic,1> y,Matrix<double,Dynamic,1> w)
//y is the ground truth, w is the weight of each parameter in the objective function
{


	return;
}
