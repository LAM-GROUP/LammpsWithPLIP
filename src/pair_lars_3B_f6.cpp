/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Aidan Thompson (SNL)
------------------------------------------------------------------------- */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pair_lars_3B_f6.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include <vector>
using namespace std;


using namespace LAMMPS_NS;

#include<iostream>

#define MAXLINE 1024
#define DELTA 4


/* ---------------------------------------------------------------------- */
double f_cut6(double rcut, double r)
{       
        double pi=3.1415926535897932;
        return 0.5*(1.0+cos(pi*(r/rcut)));
}
double f_cut6D(double rcut, double r)
{       
        double pi=3.1415926535897932;
        return -pi/rcut*0.5*sin(pi*(r/rcut));
}

double powbis6(double tmp,int l)
{
	double out=1;
	for (int i_tmp=0;i_tmp<l;i_tmp++)
		out=out*tmp;
	return out;
}

double GTO(double A,double B,double r)
{       
        if (A!=0)
                return pow(r,A)*exp(-B*r*r)/(pow(2,A)*exp(-B*r*2));
        else    
                return exp(-B*r*r)/exp(-B*2*2);



}
double GTO_D(double A,double B,double r)
{       
        if (A!=0)
                return (A-2*B*r*r)*pow(r,A-1)*exp(-B*r*r)/(pow(2,A)*exp(-B*r*2));
        else    
                return -2*B*r*exp(-B*r*r)/exp(-B*2*2);
}






PairLARS_3B_f6::PairLARS_3B_f6(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  nelements = 0;
  elements = NULL;
  nparams = maxparam = 0;
  params = NULL;
  elem2param = NULL;
  map = NULL;

  maxshort = 10;
  neighshort = NULL;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairLARS_3B_f6::~PairLARS_3B_f6()
{
  if (copymode) return;

  if (elements)
    for (int i = 0; i < nelements; i++) delete [] elements[i];
  delete [] elements;
  memory->destroy(params);
  memory->destroy(elem2param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
    delete [] map;
  }
}

/* ---------------------------------------------------------------------- */

void PairLARS_3B_f6::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

//    // two-body interactions, skip half of them
//
    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;
//
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      jtype = map[type[j]];
      ijparam = elem2param[itype][jtype][jtype];
      if (rsq >= params[ijparam].cutsq) {
        continue;
      } else {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort/2;
          memory->grow(neighshort,maxshort,"pair:neighshort");
        }
      }


      evdwl=0;
      fpair=0;


      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }

    jnumm1 = numshort - 1;

    for (jj = 0; jj < jnumm1; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      double fjxtmp,fjytmp,fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort[kk];
        ktype = map[type[k]];
        ijkparam = elem2param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        threebody(&params[ijkparam],
                  rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

        fxtmp -= fj[0] + fk[0];
        fytmp -= fj[1] + fk[1];
        fztmp -= fj[2] + fk[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLARS_3B_f6::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(neighshort,maxshort,"pair:neighshort");
  map = new int[n+1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLARS_3B_f6::settings(int narg, char **arg)
{
  if (narg != 0) error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLARS_3B_f6::coeff(int narg, char **arg)
{
  int i,j,n;

  if (!allocated) allocate();

  if (narg != 3 + atom->ntypes)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // insure I,J args are * *

  if (strcmp(arg[0],"*") != 0 || strcmp(arg[1],"*") != 0)
    error->all(FLERR,"Incorrect args for pair coefficients");

  // read args that map atom types to elements in potential file
  // map[i] = which element the Ith atom type is, -1 if NULL
  // nelements = # of unique elements
  // elements = list of element names

  if (elements) {
    for (i = 0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
  }
  elements = new char*[atom->ntypes];
  for (i = 0; i < atom->ntypes; i++) elements[i] = NULL;

  nelements = 0;
  for (i = 3; i < narg; i++) {
    if (strcmp(arg[i],"NULL") == 0) {
      map[i-2] = -1;
      continue;
    }
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    map[i-2] = j;
    if (j == nelements) {
      n = strlen(arg[i]) + 1;
      elements[j] = new char[n];
      strcpy(elements[j],arg[i]);
      nelements++;
    }
  }

  // read potential file and initialize potential parameters

  read_file(arg[2]);
  setup_params();

  // clear setflag since coeff() called once with I,J = * *

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  // set setflag i,j for type pairs where both are mapped to elements

  int count = 0;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        count++;
      }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLARS_3B_f6::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style Stillinger-Weber requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style Stillinger-Weber requires newton pair on");

  // need a full neighbor list

  int irequest = neighbor->request(this,instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLARS_3B_f6::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

void PairLARS_3B_f6::read_file(char *file)
{
  int params_per_line = 4*N_max+5;
  char **words = new char*[params_per_line+1];

  memory->sfree(params);
  params = NULL;
  nparams = maxparam = 0;

  // open file on proc 0

  FILE *fp;
  if (comm->me == 0) {
    fp = force->open_potential(file);
    if (fp == NULL) {
      char str[128];
      sprintf(str,"Cannot open Stillinger-Weber potential file %s",file);
      error->one(FLERR,str);
    }
  }

  // read each set of params from potential file
  // one set of params can span multiple lines
  // store params if all 3 element tags are in element list

  int n,nwords,ielement,jelement,kelement;
  char line[MAXLINE],*ptr;
  int eof = 0;

  while (1) {
    if (comm->me == 0) {
      ptr = fgets(line,MAXLINE,fp);
      if (ptr == NULL) {
        eof = 1;
        fclose(fp);
      } else n = strlen(line) + 1;
    }
    MPI_Bcast(&eof,1,MPI_INT,0,world);
    if (eof) break;
    MPI_Bcast(&n,1,MPI_INT,0,world);
    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // strip comment, skip line if blank

    if ((ptr = strchr(line,'#'))) *ptr = '\0';
    nwords = atom->count_words(line);
    if (nwords == 0) continue;
    // concatenate additional lines until have params_per_line words

    while (nwords < params_per_line) {
      n = strlen(line);
      if (comm->me == 0) {
        ptr = fgets(&line[n],MAXLINE-n,fp);
        if (ptr == NULL) {
          eof = 1;
          fclose(fp);
        } else n = strlen(line) + 1;
      }
      MPI_Bcast(&eof,1,MPI_INT,0,world);
      if (eof) break;
      MPI_Bcast(&n,1,MPI_INT,0,world);
      MPI_Bcast(line,n,MPI_CHAR,0,world);
      if ((ptr = strchr(line,'#'))) *ptr = '\0';
      nwords = atom->count_words(line);
    }

    cout<<nwords<<"  "<<params_per_line<<endl;
    if (nwords != params_per_line)
      error->all(FLERR,"Incorrect format in LARS Threebody potential file");

    // words = ptrs to all words in line

    nwords = 0;
    words[nwords++] = strtok(line," \t\n\r\f");
    while ((words[nwords++] = strtok(NULL," \t\n\r\f"))) continue;

    // ielement,jelement,kelement = 1st args
    // if all 3 args are in element list, then parse this line
    // else skip to next entry in file

    for (ielement = 0; ielement < nelements; ielement++)
      if (strcmp(words[0],elements[ielement]) == 0) break;
    if (ielement == nelements) continue;
    for (jelement = 0; jelement < nelements; jelement++)
      if (strcmp(words[1],elements[jelement]) == 0) break;
    if (jelement == nelements) continue;
    for (kelement = 0; kelement < nelements; kelement++)
      if (strcmp(words[2],elements[kelement]) == 0) break;
    if (kelement == nelements) continue;

    // load up parameter settings and error check their values
    if (nparams == maxparam) {
      maxparam += DELTA;
      params = (Param *) memory->srealloc(params,maxparam*(sizeof(Param)),
                                          "pair:params");
    }

    params[nparams].ielement = ielement;
    params[nparams].jelement = jelement;
    params[nparams].kelement = kelement;
    params[nparams].N_selected= atof(words[3]);
    params[nparams].cut= atof(words[4]);
    for (int i=0;i<params[nparams].N_selected;i++)
    {
    	params[nparams].epsilon_all[i]= atof(words[5+4*i]);
    	params[nparams].p_all[i]= atof(words[5+(4*i+1)]);
    	params[nparams].q_all[i]= atof(words[5+(4*i+2)]);
    	params[nparams].powerl_all[i]= atof(words[5+(4*i+3)]);
    }
    nparams++;
  }

  delete [] words;
}

/* ---------------------------------------------------------------------- */

void PairLARS_3B_f6::setup_params()
{
  int i,j,k,m,n;
  double rtmp;

  // set elem2param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem2param);
  memory->create(elem2param,nelements,nelements,nelements,"pair:elem2param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) {
            if (n >= 0) error->all(FLERR,"Potential file has duplicate entry");
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry");
        elem2param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  // (cut = a*sigma)

  for (m = 0; m < nparams; m++) {
    params[m].cutsq = params[m].cut*params[m].cut;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }
}

/* ---------------------------------------------------------------------- */


void PairLARS_3B_f6::threebody(Param *paramijk,
                       double rsq1, double rsq2,
                       double *delr1, double *delr2,
                       double *fj, double *fk, int eflag, double &eng)
{
	double fprime,gprime1,gprime2;
	double r1=sqrt(rsq1);
	double r2=sqrt(rsq2);
	double cosinus=(delr1[0]*delr2[0] + delr1[1]*delr2[1] + delr1[2]*delr2[2]) * 1.0/(r1*r2); 

	double fj_tmp[3] = { 0,0,0 };
	double fk_tmp[3] = { 0,0,0 };
	double eng_tmp=0;

	for (int i_selected=0;i_selected<paramijk->N_selected;i_selected++)
	{
		if (paramijk->epsilon_all[i_selected] !=0)		// Useless as we already know N_selected
		{
			double pow2=1;
			double pow1=1;
			if (paramijk->powerl_all[i_selected]!=0)
			{
				pow2=powbis6(cosinus,paramijk->powerl_all[i_selected]-1);
				pow1=pow2*cosinus;
			}
		
			double f1=GTO(paramijk->p_all[i_selected],paramijk->q_all[i_selected],r1 ) ;
			double f2=GTO(paramijk->p_all[i_selected],paramijk->q_all[i_selected],r2 ) ;

			double fcut1=f_cut6(paramijk->cut,r1);
			double fcut2=f_cut6(paramijk->cut,r2);

			double inv1=1.0/r1;
			double inv2=1.0/r2;


			double fr1=fcut1*f1;
			double fr2=fcut2*f2;
		
			fprime=-inv1*(GTO_D(paramijk->p_all[i_selected],paramijk->q_all[i_selected],r1)*fcut1+f1*f_cut6D(paramijk->cut,r1))*fr2*pow1;
			if (paramijk->powerl_all[i_selected]!=0)
			{
				gprime1=1.0/rsq1*     fr1*fr2*paramijk->powerl_all[i_selected]*pow1;
				gprime2=inv1*inv2*fr1*fr2*paramijk->powerl_all[i_selected]*pow2;
			}
			else
			{
				gprime1=0;
				gprime2=0;
			}
			fj_tmp[0] += paramijk->epsilon_all[i_selected]*(delr1[0]*(fprime+gprime1)-delr2[0]*gprime2);
			fj_tmp[1] += paramijk->epsilon_all[i_selected]*(delr1[1]*(fprime+gprime1)-delr2[1]*gprime2);
			fj_tmp[2] += paramijk->epsilon_all[i_selected]*(delr1[2]*(fprime+gprime1)-delr2[2]*gprime2);
		
			fprime=-inv2*( GTO_D(paramijk->p_all[i_selected],paramijk->q_all[i_selected],r2) *fcut2+f2*f_cut6D(paramijk->cut,r2))*fr1*pow1;
			if (paramijk->powerl_all[i_selected]!=0)
			{
				gprime1=1.0/rsq2*     fr2*fr1*paramijk->powerl_all[i_selected]*pow1;
				gprime2=inv1*inv2*fr2*fr1*paramijk->powerl_all[i_selected]*pow2;
			}
			else
			{
				gprime1=0;
				gprime2=0;
				
			}
			fk_tmp[0] += paramijk->epsilon_all[i_selected]*(delr2[0]*(fprime+gprime1)-delr1[0]*gprime2);
			fk_tmp[1] += paramijk->epsilon_all[i_selected]*(delr2[1]*(fprime+gprime1)-delr1[1]*gprime2);
			fk_tmp[2] += paramijk->epsilon_all[i_selected]*(delr2[2]*(fprime+gprime1)-delr1[2]*gprime2);
		
			  if (eflag) eng_tmp += paramijk->epsilon_all[i_selected]*fr1*fr2*pow1;
		}
	}
	fj[0]=fj_tmp[0];
	fj[1]=fj_tmp[1];
	fj[2]=fj_tmp[2];
	fk[0]=fk_tmp[0];
	fk[1]=fk_tmp[1];
	fk[2]=fk_tmp[2];
	if (eflag) eng=eng_tmp;

}
