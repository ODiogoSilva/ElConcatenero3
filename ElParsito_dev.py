#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ElParsito.py
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

class SeqUtils ():
	def __init__ (self, gap, missing):
		self.gap = gap
		self.missing = missing
		
	def read_alignment (self, input_alignment, alignment_format):
		""" ONLY FOR SINGLE FILE/LOCI INPUT: Function that parses an input file alignment and returns a dictionary with the taxa as keys and sequences as values """ 
		
		def rm_illegal (string):
			""" Function that removes illegal characters from taxa names """
			illegal_chars = [":",",",")","(",";","[","]","'"]
			clean_name = "".join([char for char in name if char not in illegal_chars])
			if taxa_name != clean_name:
				print ("\nWARNING: Removed illegal characters from the taxa %s" % taxa_name)
			return clean_name

		self.alignment_storage = {} # Save the taxa and their respective sequences
		self.taxa_order = [] # Save taxa names to maintain initial order
		file_handle = open(input_alignemnt)
		
		# PARSING PHYLIP FORMAT
		if alignment_format == "phylip":
			header = file_handle.readline().split() # Get the number of taxa and sequence length from the file header
			self.sequence_len = header[0],header[1]
			for line in file_handle:
				if line != "":
					taxa = line.split()[0].replace(" ","")
					taxa = rm_illegal(taxa)
					self.taxa_order.append(taxa)
					sequence = line.split()[1].strip()
					self.alignment_storage[taxa] = sequence
			
		# PARSING FASTA FORMAT
		elif alignment_format == "fasta":
			for line in file_handle:
				if line.strip().startswith(">"):
					taxa = line[1:].strip().replace(" ","_")
					taxa = rm_illegal(taxa)
					self.taxa_order.append(taxa)
					self.alignment_storage[taxa] = ""
				else:
					self.alignment_storage[taxa] += line.strip()
			self.sequence_len = len(list(self.alignment_storage.values())[0])
			
		# PARSING NEXUS FORMAT
		elif alignment_format == "nexus":
			counter = 0
			for line in fine_handle:
				if line.strip().lower() == "matrix" and counter == 0: # Skips the nexus header
					counter = 1
				elif line.strip() == ";" and counter == 1: # Stop parser here
					counter = 0
				elif line != "" and counter == 1: # Start parsing here
					taxa = line.strip().split()[0].replace(" ","")
					taxa = rm_illegal(taxa)
					self.taxa_order.append(taxa)
					if taxa in self.alignment_storage: # This accomodates for the interleave format
						self.alignment_storage[taxa] += line.strip().split()[1:]
					else:
						self.alignment_storage[taxa] = line.strip().split()[1:]
			self.sequence_len = len(list(self.alignment_storage.values())[0])
		self.check_sizes (self.alignment_storage, input_alignment)
		return (self.alignment_storage, self.taxa_order, self.sequence_len)
		
	def read_alignments (self, input_list, alignment_format):
		""" Function that parses multiple alignment/loci files and returns a dictionary with the taxa as keys and sequences as values as well as two integers corresponding to the number of taxa and sequence length """
		loci_lengths = [] # Saves the sequence lengths of the 
		self.loci_range = [] # Saves the loci names as keys and their range as values
		self.main_taxa_order = []
		for infile in input_list:
			current_alignment, taxa_order, current_sequence_len = self.read_alignment(infile,alignment_format) # Parse the current alignment
			if loci_lengths == []:
				main_alignment = current_alignment # Create the main alignment dictionary from the first current alignment and never visit this statement again
				self.main_taxa_order = taxa_order
				loci_lengths.append(current_sequence_len)
				self.loci_range.append((infile.split(".")[0],"1-%s" % (current_sequence_len-1))) # Saving the range for the first loci
			else:
				for taxa, sequence in current_alignment.items(): 
					if taxa in main_alignment: 
						main_alignment[taxa] += sequence # Append the sequence from the current alignment to the respective taxa in the main alignment
					elif taxa not in main_alignment:
						main_alignment[taxa] = "n"*sum(loci_lengths)+sequence # If the taxa does not yet exist in the main alignment, create the new entry with a sequence of 'n' characters of the same size as the length of the missed loci and the sequence from the current alignment
						self.main_taxa_order.append(taxa)
				self.loci_range.append((infile.split(".")[0],"%s-%s" % (loci_lengths[-1], current_sequence_len-1))) # Saving the range for the subsequent loci
				loci_lengths.append(current_sequence_len)
				for taxa in main_alignment:
					if taxa not in current_alignment: # Check if any taxa from the main alignment are missing from the current alignment. If yes, fill them with 'n'
						main_alignment[taxa] = "n"*current_sequence_len
		self.main_taxa_order = self.taxa_order
		return (main_alignment, self.taxa_order, self.loci_range)
		
	def check_sizes (self, alignment_dic, current_file):
		""" Function to test whether all sequences are of the same size and, if not, which are different """
		# Determine the most common length 
		commonSeq = max(set([v for v in dic.values()]),key=[v for v in dic.values()].count)
		# Creates a dictionary with the sequences, and respective length, of different length
		difLength = dict((key,value) for key, value in dic.items() if len(commonSeq) != len(value))
		if difLength != {}:
			print ("\nWARNING: Unequal sequence lenght detected in %s for the following taxa" % i_file)
			
	def zorro2rax (self, alignment_file_list, zorro_sufix="_zorro.out"):
		""" Function that converts the floating point numbers contained in the original zorro output files into intergers that can be interpreted by RAxML. If multiple alignment files are provided, it also concatenates them in the same order """
		weigths_storage = []
		for alignment_file in alignment_file_list:
			zorro_file = alignment_file.split(".")[0]+zorro_sufix # This assumes that the prefix of the alignment file is shared with the corresponding zorro file
			zorro_handle = open(zorro_file)
			weigths_storage += [round(float(weigth.strip())) for weigth in zorro_hanlde]
		return zorro_storage
	
	class writer ():
		
		# The space (in characters) available for the taxon name before the sequence begins
		self.seq_space_nex = 20
		self.seq_space_phy = 30
		self.seq_space_ima2 = 10
		# Cut the taxa names by the following character:
		self.cut_space_nex = 20
		self.cut_space_phy = 50
		self.cut_space_ima2 = 8
		
		def __init__ (self, output_file, coding):
			self.output_file = output_file
			self.coding = coding
	
		def phylip (self, alignment_dic, conversion=None):
			""" Writes a pre-parsed alignment dictionary into a new phylip file """
			out_file = open(self.output_file+".phy","w")
			out_file.write("%s %s" % (len(alignment_dic), sum(self.loci_lengths)))
			for key in self.taxa_order:
					out_file.write("%s %s\n" % (clean_key[:self.cut_space_phy].ljust(self.seq_space_phy),alignment_dic[key]))
			if conversion == None:
				partition_file = open(self.output_file+"_part.File","a")
				for partition,lrange in self.loci_range:
					partition_file.write("%s, %s = %s" % (coding,partition,lrange))
			out_file.close()
					
		def fasta (self, alignment_dic, conversion=None):
			""" Writes a pre-parsed alignment dictionary into a new fasta file """
			out_file = open(self.output_file+".fas","w")
			for key in self.taxa_order:
				out_file.write(">%s\n%s\n" % (key,alignment_dic[key]))
			out_file.close()
				
		def nexus (self, alignment_dic, conversion=None):
			""" Writes a pre-parsed alignment dictionary into a new nexus file """
			out_file = open(self.output_file+".nex","w")
			if conversion == None:
				partition_string = "".join([part[1] for part in self.loci_range])
				out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s;\n\tformat datatype=mixed (%s) interleave=no gap=%s missing=%s;\n\tmatrix\n" (len(alignment_dic), sum(self.loci_lengths), partition_string, self.gap, self.missing))
			for key in self.taxa_order:
				out_file.write("%s %s\n" % (key[:cut_space_nex].ljust(seq_space_nex),alignment_dic[key]))
			out_file.close()
				
		def zorro (self, zorro_weigths):
			""" Creates a concatenated file with the zorro weigths for the corresponding alignment files """
			outfile = self.output_file+"_zorro.out"
			for weigth in zorro_weigths:
				outfile.write("%s\n" % weigth)
			outfile.close()
