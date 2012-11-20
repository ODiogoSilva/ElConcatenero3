#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  ElParsito.py v3.1.0-0
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

def autofinder (reference_file):
	""" Autodetects the type of file to be parsed. Based on headers """
	autofind = "unknown"
	sequence = ""
	file_handle = open(reference_file,'r')
     
	header = file_handle.readline()
	while header.startswith("\n"):
		header = next(file_handle)
     
	if header.upper().startswith("#NEXUS"):
		autofind = "nexus"
		for line in file_handle:
			if line.strip().lower() == "matrix":
				sequence = "".join(file_handle.readline().split()[1:]).strip()
				break
         
	elif header.startswith(">"):
		autofind = "fasta"
		for line in file_handle:
			if line.strip()[0] != ">":
				sequence += line.strip()
			elif line.strip()[0] == ">":
				break
         
	elif len(header.strip().split()) == 2 and header[0].isdigit() and header[1].isdigit():
		autofind = "phylip"
		sequence = "".join(file_handle.readline().split()[1:]).strip()
	
	# Guessing the genetic code
	code = guess_code (sequence)
	
	return autofind, code

def guess_code (sequence):
	""" Function that guesses the code of the molecular sequences (i.e., DNA or Protein) based on the first sequence of a reference file """
	sequence = sequence.upper().replace("-","")
	DNA_count = sequence.count("A") + sequence.count("T") + sequence.count("G") + sequence.count("C") + sequence.count("N")
	DNA_proportion = float(DNA_count)/float(len(sequence))
	if DNA_proportion > 0.9:
		code = ("DNA","N")
	else:
		code = ("Protein","X")
	return code	

class SeqUtils ():
	def __init__ (self, missing="X"):
		self.missing = missing
	
	def rm_illegal (self,string):
		""" Function that removes illegal characters from taxa names """
		illegal_chars = [":",",",")","(",";","[","]","'", '"']
		clean_name = "".join([char for char in string if char not in illegal_chars])
		if string != clean_name:
			print ("\nWARNING: Removed illegal characters from the taxa %s" % string)
		return clean_name
		
	def duplicate_taxa (self, taxa_list):
		""" Function that identifies duplicated taxa """
		import collections
		duplicated_taxa = [x for x, y in collections.Counter(taxa_list).items() if y > 1]
		return duplicated_taxa
	
	def check_format (self,input_alignment,alignment_format):
		""" This function performs some very basic checks to see if the format of the input file is in accordance to the input file format specified when the script is executed """
		input_handle = open(input_alignment)
		line = input_handle.readline()
		while line.strip() == "":
			line = next(input_handle)
		
		if alignment_format == "fasta":
			if line.strip()[0] != ">":
				print ("File not in Fasta format. First non-empty line of the input file %s does not start with '>'. Please verify the file, or the input format settings\nExiting..." % input_alignment)
				raise SystemExit
		elif alignment_format == "nexus":
			if line.strip().lower() != "#nexus":
				print ("File not in Nexus format. First non-empty line of the input file %s does not start with '#NEXUS'. Please verify the file, or the input format settings\nExiting..." % input_alignment)
				raise SystemExit
		elif alignment_format == "phylip":
			try:
				header = line.strip().split()
				int(header[0])
				int(header[1])
			except:
				print ("File not in correct Phylip format. First non-empty line of the input file %s does not start with two intergers separated by whitespace. Please verify the file, or the input format settings\nExiting..." % input_alignment)
				raise SystemExit
				
	def rm_taxa (self, alignment_dic, taxa_list):
		""" Function that removes specified taxa from the alignment """
		alignment_mod = {}
		taxa_order = []
		for taxa, sequence in alignment_dic.items():
			if taxa not in taxa_list:
				alignment_mod[taxa] = sequence
				taxa_order.append(taxa)
		return alignment_mod, taxa_order
	
	def pickle_taxa (self, alignment_dic, mode):
		""" Function that exports the list of taxa from an alignment """
		import pickle
		self.taxa_list = []
		if mode == "dump":
			self.taxa_list = [taxa for taxa in alignment_dic.keys()]
			pickle.dump(self.taxa_list, open("taxa_list","wb"))
			print ("Taxa names have been saved in the pickle file 'taxa_list'\nExiting...")
			raise SystemExit
		elif mode == "load":
			self.taxa_list = pickle.load(open("taxa_list","rb"))
		return self.taxa_list
		
	def import_taxa (self, alignment_dic):
		""" Function that imports new taxa. It mainly exists to complete single locus aligments with taxa that are not present in the current alignment but occur in other alignments """
		alignment_len = self.loci_lengths[0]
		for taxa in self.taxa_list:
			if taxa not in alignment_dic:
				alignment_dic[taxa] = self.missing*alignment_len
		return alignment_dic, self.taxa_list
		
	def check_sizes (self, alignment_dic, current_file):
		""" Function to test whether all sequences are of the same size and, if not, which are different """
		# Determine the most common length 
		commonSeq = max(set([v for v in alignment_dic.values()]),key=[v for v in alignment_dic.values()].count)
		# Creates a dictionary with the sequences, and respective length, of different length
		difLength = dict((key,value) for key, value in alignment_dic.items() if len(commonSeq) != len(value))
		if difLength != {}:
			print ("\nWARNING: Unequal sequence lenght detected in %s" % current_file)

	def zorro2rax (self, alignment_file_list, zorro_sufix="_zorro.out"):
		""" Function that converts the floating point numbers contained in the original zorro output files into intergers that can be interpreted by RAxML. If multiple alignment files are provided, it also concatenates them in the same order """
		weigths_storage = []
		for alignment_file in alignment_file_list:
			zorro_file = alignment_file.split(".")[0]+zorro_sufix # This assumes that the prefix of the alignment file is shared with the corresponding zorro file
			zorro_handle = open(zorro_file)
			weigths_storage += [round(float(weigth.strip())) for weigth in zorro_handle]
		return weigths_storage

	def read_alignment (self, input_alignment, alignment_format, size_check=True):
		""" ONLY FOR SINGLE FILE/LOCI INPUT: Function that parses an input file alignment and returns a dictionary with the taxa as keys and sequences as values """ 
		
		self.check_format (input_alignment, alignment_format)

		alignment_storage = {} # Save the taxa and their respective sequences
		taxa_order = [] # Save taxa names to maintain initial order
		model = [] # Only applies for nexus format. It stores any potential substitution model at the end of the file
		
		file_handle = open(input_alignment)
		
		# PARSING PHYLIP FORMAT
		if alignment_format == "phylip":
			header = file_handle.readline().split() # Get the number of taxa and sequence length from the file header
			self.loci_lengths = int(header[1])
			for line in file_handle:
				if line != "":
					taxa = line.split()[0].replace(" ","")
					taxa = self.rm_illegal(taxa)
					taxa_order.append(taxa)
					sequence = line.split()[1].strip()
					alignment_storage[taxa] = sequence
			
		# PARSING FASTA FORMAT
		elif alignment_format == "fasta":
			for line in file_handle:
				if line.strip().startswith(">"):
					taxa = line[1:].strip().replace(" ","_")
					taxa = self.rm_illegal(taxa)
					taxa_order.append(taxa)
					alignment_storage[taxa] = ""
				else:
					alignment_storage[taxa] += line.strip()
			self.loci_lengths = len(list(alignment_storage.values())[0])
			
		# PARSING NEXUS FORMAT
		elif alignment_format == "nexus":
			counter = 0
			for line in file_handle:
				if line.strip().lower() == "matrix" and counter == 0: # Skips the nexus header
					counter = 1
				elif line.strip() == ";" and counter == 1: # Stop parser here
					counter = 2
				elif line.strip() != "" and counter == 1: # Start parsing here
					taxa = line.strip().split()[0].replace(" ","")
					taxa = self.rm_illegal(taxa)
					if taxa not in taxa_order: # Prevents duplications in the list when the input format is interleave
						taxa_order.append(taxa)
					if taxa in alignment_storage: # This accomodates for the interleave format
						alignment_storage[taxa] += "".join(line.strip().split()[1:])
					else:
						alignment_storage[taxa] = "".join(line.strip().split()[1:])
						
				# This bit of code will extract a potential substitution model from the file
				elif counter == 2 and line.lower().strip().startswith("lset"):
					model.append(line.strip())
				elif counter == 2 and line.lower().strip().startswith("prset"):
					model.append(line.strip())

			self.loci_lengths = len(list(alignment_storage.values())[0])
		
		# Checks the size consistency of the alignment
		if size_check == True:
			self.check_sizes (alignment_storage, input_alignment)
		
		# Checks for duplicate taxa
		if len(taxa_order) != len(set(taxa_order)):
			taxa = self.duplicate_taxa(taxa_order)
			print ("WARNING: Duplicated taxa have been found in file %s (%s). Please correct this problem and re-run the program\n" %(input_alignment,", ".join(taxa)))
			raise SystemExit
		
		return (alignment_storage, taxa_order, self.loci_lengths, model)
		
	def read_alignments (self, input_list, alignment_format, progress_stat=True):
		""" Function that parses multiple alignment/loci files and returns a dictionary with the taxa as keys and sequences as values as well as two integers corresponding to the number of taxa and sequence length """
		
		loci_lengths = [] # Saves the sequence lengths of the 
		loci_range = [] # Saves the loci names as keys and their range as values
		models = [] 
		main_taxa_order = []
		
		for infile in input_list:
			
			# When set to True, this statement produces a progress status on the terminal
			if progress_stat == True: 
				print ("\rProcessing file %s out of %s" % (input_list.index(infile)+1,len(input_list)),end="")

			# Parse the current alignment
			current_alignment, taxa_order, current_sequence_len, model = self.read_alignment(infile,alignment_format)
			
			# If input format is nexus, save the substution model, if any
			if alignment_format == "nexus" and model != []:
				models.append(model)

			# Algorithm that fills absent taxa with missing data
			if loci_lengths == []:
				main_alignment = current_alignment # Create the main alignment dictionary from the first current alignment and never visit this statement again
				main_taxa_order = taxa_order
				loci_lengths.append(current_sequence_len)
				loci_range.append((infile.split(".")[0],"1-%s" % (current_sequence_len))) # Saving the range for the first loci
	
			else:
				for taxa, sequence in current_alignment.items(): 
					if taxa in main_alignment: 
						main_alignment[taxa] += sequence # Append the sequence from the current alignment to the respective taxa in the main alignment
					elif taxa not in main_alignment:
						main_alignment[taxa] = self.missing*sum(loci_lengths)+sequence # If the taxa does not yet exist in the main alignment, create the new entry with a sequence of 'n' characters of the same size as the length of the missed loci and the sequence from the current alignment
						main_taxa_order.append(taxa)

				# Saving the range for the subsequent loci
				loci_range.append((infile.split(".")[0],"%s-%s" % (sum(loci_lengths)+1, sum(loci_lengths)+current_sequence_len)))
				loci_lengths.append(current_sequence_len)
				
				# Check if any taxa from the main alignment are missing from the current alignment. If yes, fill them with 'n'
				for taxa in main_alignment.keys():
					if taxa not in current_alignment: 
						main_alignment[taxa] += self.missing*current_sequence_len
		else:
			print ("\n")				
		return (main_alignment, main_taxa_order, sum(loci_lengths), loci_range, models)
	
class writer ():
		
	def __init__ (self, output_file, taxa_order, coding, loci_lengths, loci_range=None, gap = "-", missing = "n", conversion = None, models = []):
		self.output_file = output_file
		self.taxa_order = taxa_order
		self.coding = coding
		self.loci_lengths = loci_lengths
		self.loci_range = loci_range
		self.models = models
		self.gap = gap
		self.missing = missing
		# The space (in characters) available for the taxon name before the sequence begins
		self.seq_space_nex = 40
		self.seq_space_phy = 30
		self.seq_space_ima2 = 10
		# Cut the taxa names by the following character:
		self.cut_space_nex = 50
		self.cut_space_phy = 50
		self.cut_space_ima2 = 8
			
	def phylip (self, alignment_dic, conversion=None):
		""" Writes a pre-parsed alignment dictionary into a new phylip file """
		out_file = open(self.output_file+".phy","w")
		out_file.write("%s %s\n" % (len(alignment_dic), self.loci_lengths))
		for key in self.taxa_order:
				out_file.write("%s %s\n" % (key[:self.cut_space_phy].ljust(self.seq_space_phy),alignment_dic[key]))
		if conversion == None:
			partition_file = open(self.output_file+"_part.File","a")
			for partition,lrange in self.loci_range:
				partition_file.write("%s, %s = %s\n" % (self.coding,partition,lrange))
		out_file.close()
				
	def fasta (self, alignment_dic, conversion=None):
		""" Writes a pre-parsed alignment dictionary into a new fasta file """
		out_file = open(self.output_file+".fas","w")
		for key in self.taxa_order:
			out_file.write(">%s\n%s\n" % (key,alignment_dic[key]))
		out_file.close()
			
	def nexus (self, alignment_dic, conversion=None, form="leave"):
		""" Writes a pre-parsed alignment dictionary into a new nexus file """
		align_length = len(list(alignment_dic.values())[0])
		out_file = open(self.output_file+".nex","w")
		
		# This writes the output in interleave format
		if form == "interleave":
			out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=yes gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment_dic), self.loci_lengths, self.coding, self.gap, self.missing))
			counter = 0
			for i in range (90,align_length,90):
				for key in self.taxa_order:
					out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),alignment_dic[key][counter:i]))
				else:
					out_file.write("\n")
					counter = i				
			else:
				for key in self.taxa_order:
					out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),alignment_dic[key][i:align_length]))
				else:
					out_file.write("\n")
			out_file.write(";\n\tend;")
		# This writes the output in leave format (default)
		else:
			out_file.write("#NEXUS\n\nBegin data;\n\tdimensions ntax=%s nchar=%s ;\n\tformat datatype=%s interleave=no gap=%s missing=%s ;\n\tmatrix\n" % (len(alignment_dic), self.loci_lengths, self.coding, self.gap, self.missing))
			for key in self.taxa_order:
				out_file.write("%s %s\n" % (key[:self.cut_space_nex].ljust(self.seq_space_nex),alignment_dic[key]))
			out_file.write(";\n\tend;")
			
		# If this is not a conversion, write a block with the charsets and partition definition at the end of the taxa block
		if conversion == None:
			out_file.write("\nbegin mrbayes;\n")
			for partition,lrange in self.loci_range:
				out_file.write("\tcharset %s = %s;\n" % (partition,lrange))
			out_file.write("\tpartition part = %s: %s;\n\tset partition=part;\nend;\n" % (len(self.loci_range),", ".join([part[0] for part in self.loci_range])))
			
			# Concatenates the substitution models of the individual partitions
			if self.models != []:
				loci_number = 1
				out_file.write("begin mrbayes;\n")
				for model in self.models:
					m1 = model[0].split()
					m2 = model[1].split()
					m1_final = m1[0]+" applyto=("+str(loci_number)+") "+" ".join(m1[1:])
					m2_final = m2[0]+" applyto=("+str(loci_number)+") "+" ".join(m2[1:])
					out_file.write("\t%s\n\t%s\n" % (m1_final, m2_final))
					loci_number += 1
				out_file.write("end;\n")
			
		out_file.close()
			
	def zorro (self, zorro_weigths):
		""" Creates a concatenated file with the zorro weigths for the corresponding alignment files """
		outfile = self.output_file+"_zorro.out"
		outfile_handle = open(outfile,"w")
		for weigth in zorro_weigths:
			outfile_handle.write("%s\n" % weigth)
		outfile_handle.close()
