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

"""Psi4Pasture legacy psi4 modules that have been replaced within psi4-core.
Development on these modules have stopped but interested parties may install
them with this psi4 plugin wrapper.

"""


__version__ = '0.1'
__author__  = 'Psi4 Development Team'


from .locations import ccsort_sofile, transqt2_sofile
import psi4.core

psi4.core.plugin_load(ccsort_sofile)
psi4.core.plugin_load(transqt2_sofile)



