#include "Learning_mMRF.h"

Learning_mMRF::Learning_mMRF(int nlabel)
{
	_numlabel = nlabel;
	_numvar = 0;
	_para.resize(0,1);

	return;
}

int Learning_mMRF::add_cQPBF(int label, QPBpoly* A_i)
{
	cQPBFlist *temp = new cQPBFlist;
	temp->A_i = A_i;
	temp->cid = _para.rows();
	temp->next= _componentlist;
	_componentlist = temp;

	_para.resize(_para.rows()+1,1);

	if(_numvar<A_i->numvar()) _numvar = A_i->numvar();//update the number of boolean variables

	if(label<0)
	{
		for(int i=0;i<_numlabel;i++)
		{
			_lctable.addTerm2(i,_numlabel+temp->cid,1);
		}
	}
	else
	{
		_lctable.addTerm2(label,_numlabel+temp->cid,1);
	}
	
	return _componentlist->cid;
}

void Learning_mMRF::remove_cQPBF(int cid)
{
	cQPBFlist *temp = _componentlist;

	while(temp!=NULL)
	{
		if(temp->cid == cid) 
		{
				temp->A_i = NULL;
				break;
		}
	}

	for(int i=0;i<_numlabel;i++)
	{
		_lctable.delTerm2(i,_numlabel+cid);
	}

	return;
}

int Learning_mMRF::numcomp()
{
	if(_componentlist==NULL) return 0;
	return _componentlist->cid+1;
}

int Learning_mMRF::numvar()
{
	return _numvar;
}


void Learning_mMRF::learn(Matrix<int,Dynamic,1> y)
{
	Matrix<int,Dynamic,Dynamic> subid; //contains the new id of each components in each subproblem

	subid.resize(_numlabel,_para.size());

	subid.fill(0);

	Learning_QPBF *learner = new Learning_QPBF[_numlabel];

	for(cQPBFlist* templist=_componentlist;templist!=NULL;templist=templist->next)
	{
		for(int i=0;i<_numlabel;i++)
		{
			if(_lctable.getTerm2(i,_numlabel+templist->cid)>0)
			{
				subid(i,templist->cid)=learner[i].add_cQPBF(templist->A_i);
			}
		}
	}


	Matrix<double,Dynamic,Dynamic> lambda;
	lambda.resize(_numlabel,_para.size());
	lambda.fill(0);


	Matrix<double,Dynamic,1> *coeff;
	coeff = new Matrix<double,Dynamic,1>[_numlabel];

	for(int i=0;i<_numlabel;i++)
	{
		coeff[i].resize(learner[i].numcomp());

		for(int j=0;j<numcomp();j++)
		{
			if(_lctable(i,_numlabel+j)==0) continue;
			else
			{
				int count=0;
				for(int k=0;k<_numlabel;k++)
				{
					count += _lctable(k,_numlabel+j);
				}
				coeff[i](subid(i,j)) = 1.0/count;
			}
		}
	}




	return;
}

