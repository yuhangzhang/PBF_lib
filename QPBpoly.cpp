#include "QPBpoly.h"


QPBpoly::QPBpoly()
{  
	_numvar = 0;
	return;
}


void QPBpoly::addTerm1(int v0, double coeff)
{
	if(v0>=_numvar) 
	{
			_numvar = v0+1;
	}
	return addTerm2(v0,v0,coeff);
}


void QPBpoly::addTerm2(int v0, int v1, double coeff)
{
	if(v0<=v1)
	{
		if(v1>=_numvar) 
		{
			_numvar = v1+1;
		}

		std::pair<QPBF::iterator, bool> itit = _QPBcoeff.insert(std::make_pair(std::make_pair(v0,v1),coeff));

		if(itit.second==false)//already got this term
		{
			(itit.first)->second += coeff;
		}
		else 
		{
			return;
		}
	}
	else
	{
		return addTerm2(v1,v0,coeff);
	}		
}

double QPBpoly::getTerm1(int v0)
{
	return getTerm2(v0,v0);
}

double QPBpoly::getTerm2(int v0, int v1)
{
	if(v0<=v1)
	{
		QPBF::iterator itit = _QPBcoeff.find(std::make_pair(v0,v1));

		if(itit== _QPBcoeff.end())
		{
			return 0;
		}
		else
		{
			return itit->second;
		}	
	}
	else
	{
		return getTerm2(v1,v0);
	}
}



bool QPBpoly::isSubmodular()
{
	for(QPBF::iterator itit=_QPBcoeff.begin();itit!=_QPBcoeff.end();itit++)
	{
		if((itit->first).first!=(itit->first).second&&itit->second>0)
		{
			return false;
		}
	}

	return true;
}

void QPBpoly::delTerm1(int v0)
{
	return delTerm2(v0,v0);
}

void QPBpoly::delTerm2(int v0, int v1)
{
	if(v0<=v1)
	{
		_QPBcoeff.erase(std::make_pair(v0,v1));
	}
	else
	{
		delTerm2(v1,v0);
	}

	return;
}


QPBF::iterator QPBpoly::firstTerm()
{
	return _QPBcoeff.begin();
}

QPBF::iterator QPBpoly::lastTerm()
{
	return _QPBcoeff.end();
}

QPBpoly QPBpoly::operator+(QPBpoly qpbf) 
{
	int newnumvar = (_numvar>qpbf._numvar)?_numvar:qpbf._numvar;

	QPBpoly newqpbf;

	QPBF::key_compare less= _QPBcoeff.key_comp();

	QPBF::iterator it1=_QPBcoeff.begin(), it2=qpbf.firstTerm();
	
	for(;it1!=_QPBcoeff.end()&&it2!=qpbf.lastTerm();)
	{
		if(less(it1->first,it2->first))
		{
			newqpbf.addTerm2(it1->first.first,it1->first.second,it1->second);

			it1++;
		}
		else if(less(it2->first,it1->first))
		{
			newqpbf.addTerm2(it2->first.first,it2->first.second,it2->second);

			it2++;				
		}
		else
		{
			newqpbf.addTerm2(it1->first.first,it1->first.second,it1->second+it2->second);
			
			it1++;
			it2++;			
		}
	}

	for(;it1!=_QPBcoeff.end();it1++)
	{
		newqpbf.addTerm2(it1->first.first,it1->first.second,it1->second);
	}

	for(;it2!=qpbf.lastTerm();it2++)
	{
		newqpbf.addTerm2(it2->first.first,it2->first.second,it2->second);
	}



	return newqpbf;
}

QPBpoly QPBpoly::operator-(QPBpoly qpbf) 
{
	int newnumvar = (_numvar>qpbf._numvar)?_numvar:qpbf._numvar;

	QPBpoly newqpbf;

	QPBF::key_compare less= _QPBcoeff.key_comp();

	QPBF::iterator it1=_QPBcoeff.begin(), it2=qpbf.firstTerm();
	
	for(;it1!=_QPBcoeff.end()&&it2!=qpbf.lastTerm();)
	{
		if(less(it1->first,it2->first))
		{
			newqpbf.addTerm2(it1->first.first,it1->first.second,it1->second);

			it1++;
		}
		else if(less(it2->first,it1->first))
		{
			newqpbf.addTerm2(it2->first.first,it2->first.second,-(it2->second));

			it2++;				
		}
		else
		{
			if(it1->second-it2->second!=0)
			{
				newqpbf.addTerm2(it1->first.first,it1->first.second,it1->second-it2->second);
			}
			
			it1++;
			it2++;			
		}
	}

	for(;it1!=_QPBcoeff.end();it1++)
	{
		newqpbf.addTerm2(it1->first.first,it1->first.second,it1->second);
	}

	for(;it2!=qpbf.lastTerm();it2++)
	{
		newqpbf.addTerm2(it2->first.first,it2->first.second,-(it2->second));
	}

	return newqpbf;
}

QPBpoly QPBpoly::operator*(double u) 
{
	QPBpoly newqpbf;

	QPBF::key_compare less= _QPBcoeff.key_comp();

	
	for(QPBF::iterator itit=_QPBcoeff.begin();itit!=_QPBcoeff.end();itit++)
	{
		newqpbf.addTerm2(itit->first.first,itit->first.second,itit->second);
	}

	return newqpbf;
}

double& QPBpoly::operator()(int v0, int v1)
{
	int v0_,v1_;

	if(v0<=v1)
	{
		v0_ = v0;
		v1_ = v1;
	}
	else
	{
		v1_ = v0;
		v0_ = v1;
	}

	if(v1_>=_numvar) 
	{
			_numvar = v1_+1;
	}

	QPBF::iterator _Where = _QPBcoeff.lower_bound(std::make_pair(v0_,v1_));

	if (_Where == _QPBcoeff.end() || _QPBcoeff.comp(std::make_pair(v0_,v1_), _QPBcoeff._Key(_Where._Mynode())))
	{
		_Where = _QPBcoeff.insert(_Where, std::make_pair(std::make_pair(v0_,v1_), double()));
	}

	return ((*_Where).second);
}

double& QPBpoly::operator()(int v0)
{
	return operator()(v0,v0);
}


void QPBpoly::clear()
{
	return _QPBcoeff.clear();

}

void QPBpoly::clean()
{
	for(QPBF::iterator itit=_QPBcoeff.begin();itit!=_QPBcoeff.end();)
	{
		QPBF::iterator itit2=itit;
		itit++;

		if(itit2->second==0) _QPBcoeff.erase(itit2);
	}

	return;
}

int QPBpoly::size()
{
	this->clean();

	return _QPBcoeff.size();
}

double QPBpoly::evaluate(Matrix<bool,Dynamic,1> y)
{
	double value=0;

	for(QPBF::iterator itit=_QPBcoeff.begin();itit!=_QPBcoeff.end();itit++)
	{
		if(y(itit->first.first)==true&&y(itit->first.second)==true)
		{
			value+=itit->second;
		}
	}

	return value;
}

int QPBpoly::numvar()
{
	return _numvar;
}

int QPBpoly::updatenumvar(int numvar)
{
	_numvar = numvar;

	return _numvar;
}