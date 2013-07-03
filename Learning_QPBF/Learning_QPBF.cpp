#include "Learning_QPBF.h"


Learning_QPBF::Learning_QPBF()
{
	_componentlist = NULL;

	_para.resize(0,1);

	_numvar=0;


	return;
}


void Learning_QPBF::add_cQPBF(QPBpoly A_i)
{
	cQPBFlist temp;
	temp.A_i = &A_i;
	temp.cid = _para.rows();
	temp.next= _componentlist;
	_componentlist = &temp;

	_para.resize(_para.rows()+1,1);

	if(_numvar<A_i.numvar()) _numvar=A_i.numvar();//update the number of boolean variables

	return;
}



void Learning_QPBF::learn(Matrix<bool,Dynamic,1> y,Matrix<double,Dynamic,1> w)
//y is the ground truth, w is the weight of each parameter in the objective function
{
	QPBpoly vid(_numvar);

	int counter=0;

	for(int i=0;i<_numvar;i++)
	{
		vid(i,i)=counter++;//every liner term must exist, even if zero.

		for(int j=i+1;j<_numvar;j++)
		{
			for(cQPBFlist *k=_componentlist;k!=NULL;k=k->next)
			{
				if(k->A_i->getTerm2(i,j)!=0) 
				{
					vid(i,j)=counter++;
					break;
				}
			}
		}
	}


	std::vector<Triplet<double>> tripletList;
	tripletList.reserve(_numvar*_numvar*(_para.size()+4));//a rough estimation of nonzero entries in the linear constraints
	

	//start forming the contraint matrix, 
	//each row : an entry in A
	//each column : component functions, slacks, entries in P

	for(int i=0;i<_numvar;i++)
	{
		for(int j=i;j<_numvar;j++)
		{
			//EACH ROW

			if(vid(i,j)==0)
			{
			}
			else//NOT EMPTY
			{
				for(cQPBFlist *k=_componentlist;k!=NULL;k=k->next)//each component
				{
					if(k->A_i->getTerm2(i,j)!=0) tripletList.push_back(Triplet<double>(vid(i,j),k->cid,k->A_i->getTerm2(i,j)));
				}

				if(i==j) // if diagonal entry, add slack
				{
					tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+i,1));
				}

				//P_ij 
				tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+vid(i,j),-1));//P11
				tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+counter+vid(i,j),+1));//P12
				tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*2)+vid(i,j),+1));//P21
				tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*3)+vid(i,j),-1));//P22

				if(i==j) // if diagonal entry, add slack
				{
					//add 1^TP21
					for(int k=0;k<j;k++)
					{
						if(vid(k,j)!=0)
							tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*2)+vid(k,j),-2));
					}

					tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*2)+vid(i,i),-1));
					
					//add 1^TP12, which equals to the lower triangular part of P21
					
					tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*1)+vid(i,i),-1));

					for(int k=i+1;k<_numvar;k++)
					{
						if(vid(i,k)!=0)
							tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*1)+vid(i,k),-2));
					}

					//add 1^TP22
					for(int k=0;k<j;k++)
					{
						if(vid(k,j)!=0)
							tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*3)+vid(k,j),2));
					}
					for(int k=j;k<_numvar;k++)
					{
						if(vid(i,k)!=0)
							tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*3)+vid(i,k),2));
					}
				}
			}
		}
	}
	


	return;
}
