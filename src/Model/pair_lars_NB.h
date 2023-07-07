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

PairStyle(Lars_NB_xPot,PairLars_NB_xPot)

#else

#ifndef LMP_PAIR_Lars_NB_xPot_H
#define LMP_PAIR_Lars_NB_xPot_H

#include "pair_eam_power.h"

namespace LAMMPS_NS {

// need virtual public b/c of how eam_power/fs/opt inherits from it

class PairLars_NB_xPot : virtual public PairEAM_POWER {
 public:
  PairLars_NB_xPot(class LAMMPS *);
  virtual ~PairLars_NB_xPot() {}
  void coeff(int, char **);

 protected:
  void read_file(char *);
  void file2array();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: No matching element in EAM_POWER potential file

The EAM_POWER potential file does not contain elements that match the
requested elements.

E: Cannot open EAM_POWER potential file %s

The specified EAM_POWER potential file cannot be opened.  Check that the
path and name are correct.

E: Incorrect element names in EAM_POWER potential file

The element names in the EAM_POWER file do not match those requested.

*/
