
/* ---------------------------------------------------------------------- */
double fcut(double rcut, double r)
{       
        double pi=3.1415926535897932;
        return 0.5*(1.0+cos(pi*(r/rcut)));
}
double fcutD(double rcut, double r)
{       
        double pi=3.1415926535897932;
        return -pi/rcut*0.5*sin(pi*(r/rcut));
}

double powbis(double tmp,int l)
{
	double out=1;
	for (int i_tmp=0;i_tmp<l;i_tmp++)
		out=out*tmp;
	return out;
}




