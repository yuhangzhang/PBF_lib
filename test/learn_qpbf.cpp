#include "Learning_QPBF/Learning_QPBF.h"
#include "vil/vil_resample_bilin.h"
#include "vil/vil_image_view.h"
#include "vil/vil_load.h"
#include "vil/vil_save.h"



double colorDif(vil_image_view<vxl_byte>& im, int i1, int j1, int i2, int j2)
{
	double c1,c2;
	double c=0;

	for(unsigned int i=0;i<im.nplanes();i++)
	{
		c1 = im(i1,j1,i);
		c2 = im(i2,j2,i);
		c+=fabs(c1-c2);
	}

	return c/765;
}

//exe image groundtruth 
//export PATH=$PATH:/c/Program\ Files/MATLAB/R2011a/bin/win32
int main(int argc, char* argv[])
{
	vil_image_view<vxl_byte> im = vil_load(argv[1]);
	vil_image_view<vxl_byte> gt;
	vil_image_view<vxl_byte> tm = vil_load(argv[2]);
	vil_resample_bilin(tm,gt,im.ni(),im.nj());
	vil_save(gt,"o2.png");
	Learning_QPBF learner;

	QPBpoly dt_gt(im.ni()*im.nj());//component function 0: margin to ground truth
	QPBpoly dt_r(im.ni()*im.nj());//component function 1: data term on R
	QPBpoly dt_g(im.ni()*im.nj());//component function 2: data term on R
	QPBpoly dt_b(im.ni()*im.nj());//component function 3: data term on R
	QPBpoly dt_r2(im.ni()*im.nj());//component function 4: data term on R
	QPBpoly dt_g2(im.ni()*im.nj());//component function 5: data term on R
	QPBpoly dt_b2(im.ni()*im.nj());//component function 6: data term on R


	for(unsigned int i=0;i<im.ni();i++)
	{
		for(unsigned int j=0;j<im.nj();j++)
		{
			if(gt(i,j)<100)
			{
				dt_gt.addTerm1(i*im.nj()+j,-1);
			}
			else
			{
				dt_gt.addTerm1(i*im.nj()+j,1);
			}		

			dt_r.addTerm1(i*im.nj()+j,double(im(i,j,0))/255);
			dt_g.addTerm1(i*im.nj()+j,double(im(i,j,1))/255);
			dt_b.addTerm1(i*im.nj()+j,double(im(i,j,2))/255);
			dt_r2.addTerm1(i*im.nj()+j,double(im(i,j,0))*im(i,j,0)/65025);
			dt_g2.addTerm1(i*im.nj()+j,double(im(i,j,1))*im(i,j,1)/65025);
			dt_b2.addTerm1(i*im.nj()+j,double(im(i,j,2))*im(i,j,2)/65025);
		}
	}
	
	QPBpoly pt_hc(im.ni()*im.nj());//component function 7: horizontal constant edge
	QPBpoly pt_hs(im.ni()*im.nj());//component function 8: horizontal constrast sensitive edge

	for(unsigned int i=0;i<im.ni()-1;i++)
	{
		for(unsigned int j=0;j<im.nj();j++)
		{
			pt_hc.addTerm1(i*im.nj()+j,1);
			pt_hc.addTerm1(i*im.nj()+im.nj()+j,1);
			pt_hc.addTerm2(i*im.nj()+j,i*im.nj()+im.nj()+j,-2);

			double tempd =  1-colorDif(im,i,j,i+1,j);

			pt_hs.addTerm1(i*im.nj()+j,tempd);
			pt_hs.addTerm1(i*im.nj()+im.nj()+j,tempd);
			pt_hs.addTerm2(i*im.nj()+j,i*im.nj()+im.nj()+j,-2*tempd);
		}
	}

	QPBpoly pt_vc(im.ni()*im.nj());//component function 9: vertical constant edge
	QPBpoly pt_vs(im.ni()*im.nj());//component function 10: vertical constrast sensitive edge

	for(unsigned int i=0;i<im.ni();i++)
	{
		for(unsigned int j=0;j<im.nj()-1;j++)
		{
			pt_vc.addTerm1(i*im.nj()+j,1);
			pt_vc.addTerm1(i*im.nj()+j+1,1);
			pt_vc.addTerm2(i*im.nj()+j,i*im.nj()+j+1,-2);

			double tempd =  1-colorDif(im,i,j,i,j+1);

			pt_vs.addTerm1(i*im.nj()+j,tempd);
			pt_vs.addTerm1(i*im.nj()+j+1,tempd);
			pt_vs.addTerm2(i*im.nj()+j,i*im.nj()+j+1,-2*tempd);
		}
	}

	learner.add_cQPBF(&dt_gt);//dt_gt has to be the 0th one
	learner.add_cQPBF(&dt_r);
	learner.add_cQPBF(&dt_g);
	learner.add_cQPBF(&dt_b);
	learner.add_cQPBF(&dt_r2);
	learner.add_cQPBF(&dt_g2);
	learner.add_cQPBF(&dt_b2);
	learner.add_cQPBF(&pt_hc);
	learner.add_cQPBF(&pt_hs);
	learner.add_cQPBF(&pt_vc);
	learner.add_cQPBF(&pt_vs);

	Matrix<bool,Dynamic,1> y_star;

	y_star.resize(im.ni()*im.nj());

	for(unsigned int i=0;i<im.ni();i++)
	{
		for(unsigned int j=0;j<im.nj();j++)
		{
			if(gt(i,j)<100)
			{
				y_star(i*im.nj()+j) = 0;
			}
			else
			{
				y_star(i*im.nj()+j) = 1;
			}			
		}
	}

	Matrix<double,Dynamic,1> coeff;
	coeff.resize(learner.numcomp());
	coeff.fill(1);

	learner.learn(y_star,coeff);

	Matrix<bool,Dynamic,1> y;
	vil_image_view<vxl_byte> om(im.ni(),im.nj());

	learner.optimize(y);

	for(unsigned int i=0;i<om.ni();i++)
	{
		for(unsigned int j=0;j<om.nj();j++)
		{
			om(i,j) = y(i*om.nj()+j)*200;
		}
	}

	vil_save(om,"output.png");

	return 0;
}

