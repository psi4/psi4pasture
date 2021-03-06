/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/

/*! \defgroup CCSORT ccsort: Sort integrals for Coupled-Cluster Modules */

#include <string>

namespace psi { namespace ccsort {

struct Local {
  int natom;
  int nso;
  int nocc;
  int nvir;
  int domain_polar;
  int domain_mag;
  int domain_sep;
  int *aostart;
  int *aostop;
  int **domain;
  int **pairdomain;
  int *domain_len;
  int *pairdom_len;
  int *pairdom_nrlen;
  int *weak_pairs;
  double ***V;
  double ***W;
  double *eps_occ;
  double **eps_vir;
  double cutoff;
  std::string method;
  std::string weakp;
  int filter_singles;
  double weak_pair_energy;
  double cphf_cutoff;
  double core_cutoff;
  std::string freeze_core;
  std::string pairdef;
};

}} // namespace psi::ccsort