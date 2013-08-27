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

	Matrix<double,Dynamic,1> *lambda = new Matrix<double,Dynamic,1>[_numlabel];
	for(int i=0;i<_numlabel;i++)
	{
		lambda[i].resize(_learner[i].numcomp());
		lambda[i].fill(0);
	}

	double dual0=0;
	double dual1=-1;
	double bdual=0;

	while(dual0!=dual1)
	{
		dual0=dual1;
		dual1=0;

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

			dual1+=_learner[i].learn(y_i,coeff,lambda[i]);
			_learner[i].para(lambda[i]);
		}

		//update lambda
		Matrix<double,Dynamic,1> *gradient = new Matrix<double,Dynamic,1>[_numlabel];
		for(int i=0;i<_numlabel;i++)
		{
			gradient[i].resize(_learner[i].numcomp());
			gradient[i]=lambda[i];
		}

		Matrix<double,Dynamic,1> total;
		total.resize(numcomp());
		total.fill(0);

		for(QPBF::iterator qp=_lctable.firstTerm();qp!=_lctable.lastTerm();qp++)
		{
			total(int(qp->second)) +=lambda[qp->first.first](int(qp->first.second));
		}

		for(QPBF::iterator qp=_lctable.firstTerm();qp!=_lctable.lastTerm();qp++)
		{
			gradient[qp->first.first](qp->first.second) -= total(int(qp->second));
		}

		double square=0;

		for(int i=0;i<_numlabel;i++)
		{
			square+=gradient[i].squaredNorm();
		}

		if(dual1>=bdual) bdual = dual1+DELTA;
		
		double alpha = (bdual-dual1)/square;

		for(int i=0;i<_numlabel;i++)
		{
			lambda[i] += gradient[i]*alpha;
		}

	}


	Matrix<double,Dynamic,1> count;
	count.resize(_para.size());
	count.fill(0);
	for(QPBF::iterator qp=_lctable.firstTerm();qp!=_lctable.lastTerm();qp++)
	{
		_para(int(qp->second))+=lambda[qp->first.first](int(qp->first.second));
		count(int(qp->second)) += 1;
	}

	for(int i=0;i<_para.size();i++)
	{
		_para(i) /= count(i);
	}

	return;
}

