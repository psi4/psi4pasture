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
/*
**  CCSORT: Program to reorganize integrals for CC and MBPT calculations.
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#include "globals.h"

namespace psi { namespace ccsort {

#define IOFF_MAX 32641

void init_io();
void init_ioff(void);
void title(void);
void get_params(Options & options);
void get_moinfo(void);
void sort_oei(void);
void sort_tei(void);
void b_sort(void);
void c_sort(void);
void d_sort(void);
void e_sort(void);
void f_sort(void);
void d_spinad(void);
void e_spinad(void);
void f_spinad(void);
void scf_check(void);
void fock(void);
void denom(void);
void exit_io(void);
void cleanup(void);
int **cacheprep_uhf(int level, int *cachefiles);
int **cacheprep_rhf(int level, int *cachefiles);
void cachedone_uhf(int **cachelist);
void cachedone_rhf(int **cachelist);
void local_init(Options & options);
void local_done(void);
void cc_memcheck(void);

extern "C"
int read_options(std::string name, Options& options)
{
  if(name == "CCSORT"|| options.read_globals()) {
      /*- MODULEDESCRIPTION Sorts integrals for efficiency. Called before (non density-fitted) MP2 and
          coupled cluster computations. -*/
    // Temporary: turn on/off cctransort module.  Remove after stability is proven. -TDC (1/19/2016)
    options.add_bool("RUN_CCTRANSORT", true);
    /*- Wavefunction type !expert -*/
    options.add_str("WFN", "");
    /*- Reference wavefunction type -*/
    options.add_str("REFERENCE", "RHF");
    /*- Reference wavefunction type for EOM computations -*/
    options.add_str("EOM_REFERENCE","RHF");
    /*- The response property desired.  The unique acceptable values is ``POLARIZABILITY``
    for dipole-polarizabilitie. -*/
    options.add_str("PROPERTY", "POLARIZABILITY");
    /*- Do simulate the effects of local correlation techniques? -*/
    options.add_bool("LOCAL", false);
    /*- Value (always between one and zero) for the Broughton-Pulay completeness
    check used to contruct orbital domains for local-CC calculations. See
    J. Broughton and P. Pulay, J. Comp. Chem. 14, 736-740 (1993) and C. Hampel
    and H.-J. Werner, J. Chem. Phys. 104, 6286-6297 (1996). -*/
    options.add_double("LOCAL_CUTOFF", 0.02);
    /*- Cutoff value for local-coupled-perturbed-Hartree-Fock -*/
    options.add_double("LOCAL_CPHF_CUTOFF", 0.10);
    /*- Local core cutoff value -*/
    options.add_double("LOCAL_CORE_CUTOFF",0.05);
    /*- Type of local-CCSD scheme to be simulated. ``WERNER`` (unique avaliable option) selects the method
    developed by H.-J. Werner and co-workers. -*/
    options.add_str("LOCAL_METHOD","WERNER");
    /*- Desired treatment of "weak pairs" in the local-CCSD method. The value of ``NONE`` (unique avaliable option) treats weak pairs in
    the same manner as strong pairs. -*/
    options.add_str("LOCAL_WEAKP","NONE");
    /*- Definition of local pair domains, unique avaliable option is BP, Boughton-Pulay. -*/
    options.add_str("LOCAL_PAIRDEF","BP");
    /*- Do use augment domains with polarized orbitals? -*/
    options.add_bool("LOCAL_DOMAIN_POLAR", false);
    /*- Do generate magnetic-field CPHF solutions for local-CC? -*/
    options.add_bool("LOCAL_DOMAIN_MAG", false);
    /*- -*/
    options.add_bool("LOCAL_DOMAIN_SEP", false);
    /*- Do apply local filtering to single excitation amplitudes? -*/
    options.add_bool("LOCAL_FILTER_SINGLES", false);
    /*- The algorithm to use for the $\left<VV||VV\right>$ terms -*/
    options.add_str("AO_BASIS", "NONE", "NONE DISK DIRECT");
    /*- Do retain the input two-electron integrals? -*/
    options.add_bool("KEEP_TEIFILE", false);
    /*- Do retain the input one-electron integrals? -*/
    options.add_bool("KEEP_OEIFILE", false);
    /*- Minimum absolute value below which integrals are neglected. -*/
    options.add_double("INTS_TOLERANCE", 1e-14);
    /*- Cacheing level for libdpd governing the storage of amplitudes,
    integrals, and intermediates in the CC procedure. A value of 0 retains
    no quantities in cache, while a level of 6 attempts to store all
    quantities in cache.  For particularly large calculations, a value of
    0 may help with certain types of memory problems.  The default is 2,
    which means that all four-index quantites with up to two virtual-orbital
    indices (e.g., $\langle ij | ab \rangle>$ integrals) may be held in the cache. -*/
    options.add_int("CACHELEVEL", 2);
    /*- Energy of applied field [au] for dynamic properties -*/
    options.add("OMEGA", new ArrayType());
    /*- Convert ROHF MOs to semicanonical MOs -*/
    options.add_bool("SEMICANONICAL", true);
  }
  return true;
}

extern "C"
SharedWavefunction ccsort(SharedWavefunction ref_wfn, Options &options)
{

  int i;
  int **cachelist, *cachefiles;
  int bamount,famount; /* Truncated theoretical number of B/F-type ints that
                          could be stored in cache at once                  */

  unsigned long int ia_size, ab_size, ij_size, f_size, t2_size, b_size;

  init_io();
  init_ioff();
  title();

  get_params(options);
  get_moinfo();

  cachefiles = init_int_array(PSIO_MAXUNIT);

  if(params.ref == 2) { /*** UHF references ***/
    cachelist = cacheprep_uhf(params.cachelev, cachefiles);

    std::vector<int*> spaces;
    spaces.push_back(moinfo.aoccpi);
    spaces.push_back(moinfo.aocc_sym);
    spaces.push_back(moinfo.avirtpi);
    spaces.push_back(moinfo.avir_sym);
    spaces.push_back(moinfo.boccpi);
    spaces.push_back(moinfo.bocc_sym);
    spaces.push_back(moinfo.bvirtpi);
    spaces.push_back(moinfo.bvir_sym);
    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 4, spaces);
  }
  else { /*** RHF/ROHF references ***/
    cachelist = cacheprep_rhf(params.cachelev, cachefiles);

    std::vector<int*> spaces;
    spaces.push_back(moinfo.occpi);
    spaces.push_back(moinfo.occ_sym);
    spaces.push_back(moinfo.virtpi);
    spaces.push_back(moinfo.vir_sym);
    dpd_init(0, moinfo.nirreps, params.memory, 0, cachefiles, cachelist, NULL, 2, spaces);
  }

  /* run a small computation of memory and disk requirements */
  cc_memcheck();

  outfile->Printf( "\n");

  sort_oei();
  sort_tei();
  c_sort();
  d_sort();
  e_sort();
  f_sort();
  if(params.ref == 0) {
    d_spinad();
    e_spinad();
/*     f_spinad(); */
  }
  scf_check();
  fock();
  denom();

  /* CPHF stuff for local correlation tests */
  if(params.local) {
    local_init(options);
    local_done();
  }

  dpd_close(0);

  if(params.ref == 2) cachedone_uhf(cachelist);
  else cachedone_rhf(cachelist);
  free(cachefiles);

  cleanup();
  exit_io();
  return (ref_wfn);
}

void init_io()
{
//  int i, num_unparsed;
//  char **argv_unparsed;

//  argv_unparsed = (char **) malloc(argc * sizeof(char *));
  params.reset = 0;
//  for(i=1, num_unparsed=0; i < argc; i++) {
//    if(!strcmp(argv[i], "--reset")) params.reset = 1;
//    else argv_unparsed[num_unparsed++] = argv[i];
//  }

//  free(argv_unparsed);

  tstart();

  int i;
  if(params.reset) for(i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i,0);
  else for(i=PSIF_CC_MIN; i <= PSIF_CC_MAX; i++) psio_open(i,1);
}

void title(void)
{
  outfile->Printf( "\n");
  outfile->Printf( "\t\t\t**************************\n");
  outfile->Printf( "\t\t\t*                        *\n");
  outfile->Printf( "\t\t\t*         CCSORT         *\n");
  outfile->Printf( "\t\t\t*                        *\n");
  outfile->Printf( "\t\t\t**************************\n");
  outfile->Printf( "\n");
}

void exit_io(void)
{
  int i;
  for(i=PSIF_CC_MIN; i < PSIF_CC_TMP; i++) psio_close(i,1);
  for(i=PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) psio_close(i,0);  /* get rid of TMP files */
  for(i=PSIF_CC_TMP11+1; i <= PSIF_CC_MAX; i++) psio_close(i,1);

  tstop();
}

void init_ioff(void)
{
  int i;
  ioff = init_int_array(IOFF_MAX);
  ioff[0] = 0;
  for(i=1; i < IOFF_MAX; i++) ioff[i] = ioff[i-1] + i;
}

void cleanup(void)
{
  int i;

  psio_write_entry(PSIF_CC_INFO, "Reference Energy", (char *) &(moinfo.eref),
           sizeof(double));

  if(params.ref == 2) {

    free(moinfo.pitz2qt_A);
    free(moinfo.pitz2qt_B);
    free(moinfo.qt2pitz_A);
    free(moinfo.qt2pitz_B);

    free(moinfo.aocc);
    free(moinfo.bocc);
    free(moinfo.avir);
    free(moinfo.bvir);
    free(moinfo.all_aocc);
    free(moinfo.all_bocc);
    free(moinfo.all_avir);
    free(moinfo.all_bvir);
    free(moinfo.aoccpi);
    free(moinfo.boccpi);
    free(moinfo.avirtpi);
    free(moinfo.bvirtpi);
    free(moinfo.all_aoccpi);
    free(moinfo.all_boccpi);
    free(moinfo.all_avirtpi);
    free(moinfo.all_bvirtpi);

    free(moinfo.cc_aocc);
    free(moinfo.cc_bocc);
    free(moinfo.cc_avir);
    free(moinfo.cc_bvir);
    free(moinfo.qt_aocc);
    free(moinfo.qt_bocc);
    free(moinfo.qt_avir);
    free(moinfo.qt_bvir);
    free(moinfo.aocc_sym);
    free(moinfo.bocc_sym);
    free(moinfo.avir_sym);
    free(moinfo.bvir_sym);

    free(moinfo.cc_allaocc);
    free(moinfo.cc_allbocc);
    free(moinfo.cc_allavir);
    free(moinfo.cc_allbvir);
    free(moinfo.qt_allaocc);
    free(moinfo.qt_allbocc);
    free(moinfo.qt_allavir);
    free(moinfo.qt_allbvir);
    free(moinfo.allaocc_sym);
    free(moinfo.allbocc_sym);
    free(moinfo.allavir_sym);
    free(moinfo.allbvir_sym);

    free(moinfo.aocc_off);
    free(moinfo.bocc_off);
    free(moinfo.avir_off);
    free(moinfo.bvir_off);
    free(moinfo.all_aocc_off);
    free(moinfo.all_bocc_off);
    free(moinfo.all_avir_off);
    free(moinfo.all_bvir_off);
  }
  else {

    free(moinfo.pitz2qt);
    free(moinfo.qt2pitz);

    free(moinfo.occ);
    free(moinfo.vir);
    free(moinfo.all_occ);
    free(moinfo.all_vir);
    free(moinfo.socc);
    free(moinfo.all_socc);
    free(moinfo.occpi);
    free(moinfo.virtpi);
    free(moinfo.all_occpi);
    free(moinfo.all_virtpi);
    free(moinfo.orbsym);

    free(moinfo.cc_occ);
    free(moinfo.cc_vir);
    free(moinfo.qt_occ);
    free(moinfo.qt_vir);
    free(moinfo.occ_sym);
    free(moinfo.vir_sym);

    free(moinfo.cc_allocc);
    free(moinfo.cc_allvir);
    free(moinfo.qt_allocc);
    free(moinfo.qt_allvir);
    free(moinfo.allocc_sym);
    free(moinfo.allvir_sym);

    free(moinfo.occ_off);
    free(moinfo.vir_off);
    free(moinfo.all_occ_off);
    free(moinfo.all_vir_off);

  }

  free(moinfo.pitzer2qt);
  free(moinfo.qt2pitzer);

  // Wavefunction owns these arrays
//  free(moinfo.sopi);
//  free(moinfo.orbspi);
  free(moinfo.clsdpi);
//  free(moinfo.openpi);
//  free(moinfo.uoccpi);
//  free(moinfo.fruocc);
//  free(moinfo.frdocc);
  for(i=0; i < moinfo.nirreps; i++)
    free(moinfo.labels[i]);
  free(moinfo.labels);
  free(moinfo.frozen);
  free(ioff);
}

}} //namespace psi::ccsort
