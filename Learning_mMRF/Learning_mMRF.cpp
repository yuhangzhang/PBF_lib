#include "Learning_mMRF.h"

Learning_mMRF::Learning_mMRF(int nlabel)
{
	_numlabel = nlabel;
	_numvar = 0;
	_para.resize(0,1);

	_coeff.fill(0);

	_learner = new Learning_QPBF[_numlabel];

	return;
}

int Learning_mMRF::add_cQPBF(int label, QPBpoly* A_i)
{
	if(_numvar<A_i->numvar()) _numvar = A_i->numvar();//update the number of boolean variables

	if(label<0)
	{
		for(int i=0;i<_numlabel;i++)
		{
			_lctable.addTerm2(i,_numlabel+_learner[i].add_cQPBF(A_i),_para.size());
		}
		_coeff(_para.size())=1.0/_numlabel;
	}
	else
	{
		_lctable.addTerm2(label,_numlabel+_learner[label].add_cQPBF(A_i),_para.size());
		_coeff(_para.size())=1.0;
	}
	
	_para.resize(_para.size()+1);

	return _para.size()-1;
}

//void Learning_mMRF::remove_cQPBF(int cid)
//{
//	cQPBFlist *temp = _componentlist;
//
//	while(temp!=NULL)
//	{
//		if(temp->cid == cid) 
//		{
//				temp->A_i = NULL;
//				break;
//		}
//	}
//
//	for(int i=0;i<_numlabel;i++)
//	{
//		_lctable.delTerm2(i,_numlabel+cid);
//	}
//
//	return;
//}

int Learning_mMRF::numcomp()
{
	return _para.size();
}

int Learning_mMRF::numvar()
{
	return _numvar;
}


void Learning_mMRF::learn(Matrix<int,Dynamic,1> y)
{

	Matrix<double,Dynamic,Dynamic> lambda;
	lambda.resize(_numlabel,_para.size());
	lambda.fill(0);

	for(int i=0;i<_numlabel;i++)
	{
		Matrix<double,Dynamic,1> coeff;
		coeff.resize(_learner[i].numcomp());

		Matrix<bool,Dynamic,1> y_i;

		for(int j=0;j<y.size();j++)
		{
			y_i(j) = (y(j)==i)?1:0;
		}

		for(int j=0;j<coeff.size();j++)
		{
			coeff(j) = _coeff(int(_lctable.getTerm2(i,_numlabel+j)));
		}

		_learner[i].learn(y_i,coeff,lambda.row(i));
	}






	return;
}

