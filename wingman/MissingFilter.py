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

	def __init__ (self, alignment_object, gap_threshold=50, missing_threshold=75, gap_symbol="-", missing_symbol="n"):
		""" the gap_threshold variable is a cut-off to total_missing_proportion and missing_threshold in a cut-off to missing_proportion """

		self.alignment_obj = alignment_object
		self.alignment = alignment_object.alignment
		self.gap = gap_symbol
		self.missing = missing_symbol

		# Definig thresholds
		self.gap_threshold = gap_threshold
		self.missing_threshold = missing_threshold

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

	def filter_columns (self):
		""" Here several missing data metrics are calculated, and based on some user defined thresholds, columns with inappropriate missing data are removed """

		def delete_column (column_position):
			""" Converts alignment strings into lists and removes the ith column from the alignment """ 

			for taxa,seq in self.alignment.items():

				new_seq = list(seq)
				del(new_seq[column_position])

				self.alignment[taxa] = "".join(new_seq)

			return 0

		taxa_number = len(self.alignment)

		# Creating the column list variable
		for column_position in range(self.alignment_obj.locus_length-1, -1, -1): # The reverse iteration over the sequences is necessary to maintain the column numbers when removing them

			column = [char[column_position] for char in self.alignment]

			# Calculating metrics
			gap_proportion = (float(column.count(self.gap))/float(taxa_number))*float(100)
			missing_proportion = (float(column.count(self.missing))/float(taxa_number))*float(100)
			total_missing_proportion = gap_proportion+missing_proportion

			if total_missing_proportion > float(self.gap_threshold):

				delete_column (column_position)

			elif missing_proportion > float(self.missing_threshold):

				delete_column (column_position)