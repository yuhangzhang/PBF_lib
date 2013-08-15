#include "Learning_QPBF.h"


Learning_QPBF::Learning_QPBF()
{
	_componentlist = NULL;

	_para.resize(0,1);

	_numvar=0;


	return;
}


void Learning_QPBF::add_cQPBF(QPBpoly* A_i)
{
	cQPBFlist *temp = new cQPBFlist;
	temp->A_i = A_i;
	temp->cid = _para.rows();
	temp->next= _componentlist;
	_componentlist = temp;

	_para.resize(_para.rows()+1,1);

	//if(_numvar>0) 
	//{
	//	printf("_numvar=%d %f %f\n",_numvar,_componentlist->A_i->getTerm1(0),((_componentlist->next)->A_i)->getTerm1(0));
	//	getchar();
	//}

	if(_numvar<A_i->numvar()) _numvar=A_i->numvar();//update the number of boolean variables

	

	return;
}



void Learning_QPBF::learn(Matrix<bool,Dynamic,1> y)
//y is the ground truth, w is the weight of each parameter in the objective function
{
	QPBpoly vid(_numvar);

	int counter=1;//we count from 1 because will use vid(i,j)=0 to identify absent terms
				  //also, this simplifies the parameter transition to Matlab

	for(cQPBFlist *k=_componentlist;k!=NULL;k=k->next)
	{
		for(QPBF::iterator it=k->A_i->firstTerm();it!=k->A_i->lastTerm();it++)
		{
			if(vid.getTerm2(it->first.first,it->first.second)==0)
			{
				vid(it->first.first,it->first.second)=counter++;
			}
		}
	}

	//counter = number of active terms+1


	std::vector<Triplet<double>> tripletList;
	tripletList.reserve(4*_numvar*(_para.size()+4));//a rough estimation of nonzero entries in the linear constraints
	
	printf("Construct Matrix A ...\n");
	//start forming the contraint matrix, 
	//each row : an entry in A
	//each column : component functions, slacks, entries in P
	//however, do not do it entry by entry. Use the sparsity
	for(cQPBFlist *k=_componentlist;k!=NULL;k=k->next)
	{
		for(QPBF::iterator it=k->A_i->firstTerm();it!=k->A_i->lastTerm();it++)
		{
			int i=it->first.first;
			int j=it->first.second;
			
			tripletList.push_back(Triplet<double>(vid(i,j),k->cid+1,it->second));//remember we count from 1

		}
	}

	for(QPBF::iterator it=vid.firstTerm();it!=vid.lastTerm();it++)
	{		
		int i=it->first.first;
		int j=it->first.second;

		
		if(i==j) // if diagonal entry
		{
			// add slack
			if(y(i)==0) tripletList.push_back(Triplet<double>(it->second,_para.rows()+i+1,1));
			else tripletList.push_back(Triplet<double>(it->second,_para.rows()+i+1,-1));
			//add 1^TP21
			tripletList.push_back(Triplet<double>(it->second,_para.rows()+_numvar+(counter*2)+it->second+1,-1));
			//add 1^TP12, which equals to the lower triangular part of P21
			tripletList.push_back(Triplet<double>(it->second,_para.rows()+_numvar+(counter*1)+it->second+1,-1));
		}	

		
		if(i<j)
		{
			//add 1^TP21
			tripletList.push_back(Triplet<double>(vid.getTerm1(j),_para.rows()+_numvar+(counter*2)+it->second+1,-2));
			//add 1^TP22
			tripletList.push_back(Triplet<double>(vid.getTerm1(j),_para.rows()+_numvar+(counter*3)+it->second+1,2));

		}
		
		if(i>j)
		{
			//add 1^TP12, which equals to the lower triangular part of P21
			tripletList.push_back(Triplet<double>(vid.getTerm1(i),_para.rows()+_numvar+(counter*1)+it->second+1,-2));
			//add 1^TP22
			tripletList.push_back(Triplet<double>(vid.getTerm1(i),_para.rows()+_numvar+(counter*3)+it->second+1,2));

		}

		tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+it->second+1,-1));//P11
		tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+counter+it->second+1,+1));//P12
		tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*2)+it->second+1,+1));//P21
		tripletList.push_back(Triplet<double>(vid(i,j),_para.rows()+_numvar+(counter*3)+it->second+1,-1));//P22
	}



	printf("make f_p(y)=0\n");
	for(QPBF::iterator it=vid.firstTerm();it!=vid.lastTerm();it++)
	{
		int i=it->first.first;
		int j=it->first.second;

		if(y(i)==0&&y(j)==0)
		{
			tripletList.push_back(Triplet<double>(counter,_para.rows()+_numvar+it->second+1,1));//P11
		}
		else if(y(i)==0)
		{
			tripletList.push_back(Triplet<double>(counter,_para.rows()+_numvar+(counter*2)+it->second+1,1));//P21
		}
		else if(y(j)==0)
		{
			tripletList.push_back(Triplet<double>(counter,_para.rows()+_numvar+(counter*1)+it->second+1,1));//P12
		}
		else
		{
			tripletList.push_back(Triplet<double>(counter,_para.rows()+_numvar+(counter*3)+it->second+1,1));//P22
		}
	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Call Matlab Engine
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	printf("passing to matlab\n");
	Engine *eg = engOpenSingleUse(NULL, NULL, NULL);
	engEvalString(eg,"cd C:\\yuhang\\coin-or\\QPBF_learning\\Release\\cvx\\cvx;");
	engEvalString(eg,"cvx_setup;");	
	engEvalString(eg,"clear all;");	

	char *screen = new char[1000000];
	engOutputBuffer(eg, screen, 1000000);


	mxArray *scalar = mxCreateDoubleMatrix(1,1,mxREAL);

	*((double *) mxGetPr(scalar))=tripletList.size();
	engPutVariable(eg,"nznode",scalar);
	*((double *) mxGetPr(scalar))=counter;
	engPutVariable(eg,"numact",scalar);
	*((double *) mxGetPr(scalar))=_para.size();
	engPutVariable(eg,"numcom",scalar);
	*((double *) mxGetPr(scalar))=_numvar;
	engPutVariable(eg,"numvar",scalar);

	mxArray *index_i = mxCreateDoubleMatrix(tripletList.size(),1,mxREAL);
	mxArray *index_j = mxCreateDoubleMatrix(tripletList.size(),1,mxREAL);
	mxArray *value_s = mxCreateDoubleMatrix(tripletList.size(),1,mxREAL);

	for(int i=0;i<tripletList.size();i++)
	{
		mxGetPr(index_i)[i] = tripletList[i].row();
		mxGetPr(index_j)[i] = tripletList[i].col();
		mxGetPr(value_s)[i] = tripletList[i].value();
	}

	tripletList.clear();

	engPutVariable(eg,"index_i",index_i);
	//engEvalString(eg,"index_i = index_i+1;");
	engPutVariable(eg,"index_j",index_j);
	//engEvalString(eg,"index_j = index_j+1;");
	engPutVariable(eg,"value_s",value_s);

	engEvalString(eg,"A = sparse(index_i,index_j,value_s,numact+1,numcom+numvar+numact*4,nznode);");
	
	engEvalString(eg,"cvx_begin");

	engEvalString(eg,"variable d(numvar);");
	engEvalString(eg,"variable w(numcom);");
	engEvalString(eg,"variable p(numact*4);");


	engEvalString(eg,"minimize(norm(w,1)+norm(d,1));");
	engEvalString(eg,"subject to");
	engEvalString(eg,"A*[w;d;p]==zeros(numact+1,1);");
	engEvalString(eg,"w(0)==1;");
	engEvalString(eg,"d>=zeros(numvar,1);");
	engEvalString(eg,"p>=zeros(numact*4,1);");
	engEvalString(eg,"cvx_end");
	
	mxArray *slacks = engGetVariable(eg,"d");
	mxArray *weights = engGetVariable(eg,"w");
	mxArray *positerms = engGetVariable(eg,"p");

	
	for(int i=0;i<_para.size();i++)
	{
		_para(i) = mxGetPr(weights)[i];
	}

	return;
}
