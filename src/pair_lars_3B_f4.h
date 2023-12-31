/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
  DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lars_3B_f4,PairLARS_3B_f4)

#else

#ifndef LMP_PAIR_LARS_3B_f4_H
#define LMP_PAIR_LARS_3B_f4_H

#include "pair.h"
#include <vector>
using namespace std;

namespace LAMMPS_NS {

class PairLARS_3B_f4 : public Pair {
 public:
  PairLARS_3B_f4(class LAMMPS *);
  virtual ~PairLARS_3B_f4();
  virtual void compute(int, int);
  void settings(int, char **);
  virtual void coeff(int, char **);
  virtual double init_one(int, int);
  virtual void init_style();

  struct Param {
    int N_selected;
    double epsilon_all[100];
    double p_all[100];
    double q_all[100];
    double powerl_all[100];
    double cut,cutsq;
    int ielement,jelement,kelement;
    double epsilon,p,q;
    int powerl;

//    double omega,p,q;
//    int l;
  };

 protected:
  int N_max=35; 
  double cutmax;                // max cutoff for all elements
  int nelements;                // # of unique elements
  char **elements;              // names of unique elements
  int ***elem2param;            // mapping from element triplets to parameters
  int *map;                     // mapping from atom types to elements
  int nparams;                  // # of stored parameter sets
  int maxparam;                 // max # of parameter sets
  Param *params;                // parameter set for an I-J-K interaction
  int maxshort;                 // size of short neighbor list array
  int *neighshort;              // short neighbor list array

  virtual void allocate();
  void read_file(char *);
  virtual void setup_params();
  void threebody(Param *, double, double, double *, double *,
                 double *, double *, int, double &);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style Stillinger-Weber requires atom IDs

This is a requirement to use the LARS_3B_f4 potential.

E: Pair style Stillinger-Weber requires newton pair on

See the newton command.  This is a restriction to use the LARS_3B_f4
potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

E: Cannot open Stillinger-Weber potential file %s

The specified LARS_3B_f4 potential file cannot be opened.  Check that the path
and name are correct.

E: Incorrect format in Stillinger-Weber potential file

Incorrect number of words per line in the potential file.

E: Illegal Stillinger-Weber parameter

One or more of the coefficients defined in the potential file is
invalid.

E: Potential file has duplicate entry

The potential file has more than one entry for the same element.

E: Potential file is missing an entry

The potential file does not have a needed entry.

*/
