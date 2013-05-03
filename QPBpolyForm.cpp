#include "QPBpolyForm.h"


QPBpolyForm::QPBpolyForm()
{  
	return;
}


void QPBpolyForm::addTerm1(int v0, double coeff)
{
	return addTerm2(v0,v0,coeff);
}


void QPBpolyForm::addTerm2(int v0, int v1, double coeff)
{
	if(v1>_maxvar) _maxvar = v1;

	if(v0<=v1)
	{
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

double QPBpolyForm::getTerm2(int v0, int v1)
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



bool QPBpolyForm::isSubmodular()
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

void QPBpolyForm::delTerm1(int v0)
{
	return delTerm2(v0,v0);
}

void QPBpolyForm::delTerm2(int v0, int v1)
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


QPBF::iterator QPBpolyForm::firstTerm()
{
	return _QPBcoeff.begin();
}

QPBF::iterator QPBpolyForm::lastTerm()
{
	return _QPBcoeff.end();
}

QPBpolyForm QPBpolyForm::operator+(QPBpolyForm qpbf) 
{
	QPBpolyForm newqpbf;

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

QPBpolyForm QPBpolyForm::operator-(QPBpolyForm qpbf) 
{
	QPBpolyForm newqpbf;

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


double& QPBpolyForm::operator()(int v0, int v1)
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

	QPBF::iterator _Where = _QPBcoeff.lower_bound(std::make_pair(v0_,v1_));

	if (_Where == _QPBcoeff.end() || _QPBcoeff.comp(std::make_pair(v0_,v1_), _QPBcoeff._Key(_Where._Mynode())))
	{
		_Where = _QPBcoeff.insert(_Where, std::make_pair(std::make_pair(v0_,v1_), double()));
	}

	return ((*_Where).second);
}

double& QPBpolyForm::operator()(int v0)
{
	return operator()(v0,v0);
}


void QPBpolyForm::clear()
{
	_QPBcoeff.clear();
}

void QPBpolyForm::clean()
{
	for(QPBF::iterator itit=_QPBcoeff.begin();itit!=_QPBcoeff.end();)
	{
		QPBF::iterator itit2=itit;
		itit++;

		if(itit2->second==0) _QPBcoeff.erase(itit2);
	}

	return;
}

int QPBpolyForm::size()
{
	this->clean();

	return _QPBcoeff.size();
}
