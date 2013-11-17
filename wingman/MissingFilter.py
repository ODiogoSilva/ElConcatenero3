#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  
#  Copyright 2012 Unknown <diogo@arch>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

class MissingFilter ():
	""" Contains several methods used to trim and filter missing data from alignments. It's mainly used for inheritance """

	def __init__ (self, alignment_dict, gap_symbol="-", missing_symbol="n"):

		self.alignment = alignment_dict
		self.gap = gap_symbol
		self.missing = missing_symbol

	def filter_terminals (self):
		""" Given an alignment, this will replace the gaps in the extremities of the alignment with missing data """

		for taxa,seq in self.alignment.items():

			trim_seq = list(seq)
			counter, reverse_counter = 0, -1

			while trim_seq[counter] == self.gap:
				trim_seq[counter] = self.missing

			while trim_seq[reverse_counter] == self.gap:
				trim_seq[reverse_counter] = self.missing

			seq = "".join(trim_seq)

			self.alignment[taxa] = seq