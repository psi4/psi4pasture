#
# @BEGIN LICENSE
#
# Psi4: an open-source quantum chemistry software package
#
# Copyright (c) 2007-2016 The Psi4 Developers.
#
# The copyrights for code used from other parties are included in
# the corresponding files.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# @END LICENSE
#

"""
This module contains functions that define runtime module execution sequences
for legacy modules within psi4pasture.
"""

import psi4.core
from .locations import ccsort_sofile,transqt2_sofile

def ccsort_transqt2(ref_wfn):
    """ transqt2/ccsort runtime module execution sequence (replacement for cctransort)

    :examples:
    >>> # [1] Coupled-cluster singles and doubles calculation
    >>> # using transqt2/ccsort inplace of cctransort

    >>> molecule H2 {\n0 1\nH \nH 1 0.74\n}
    >>> set basis cc-PVDZ
    >>> set run_cctransort false
    >>> energy('ccsd')
    """
    psi4.core.plugin(transqt2_sofile,ref_wfn)
    psi4.core.plugin(ccsort_sofile,ref_wfn)
