#!/usr/bin/python3

# ElConcatenero v3.1.6
# Author: Diogo N Silva
# Last update: 10/08/2013
# ElConcatenero is tool to convert and concatenate several commonly used data format types. Currently, supported input formats include Nexus, FastA and Phylip. Output may be in Nexus, Phylip (wiht part file for RaXML) or FastA. Please type "ElConcatenero -h" or read the README.md file for information on usage.

#  Copyright 2012 Diogo N Silva <diogo@arch>
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
#  along with this program; if not,  If not, see <http://www.gnu.org/licenses/>.

import argparse
#import ElParsito3 as ep
from wingman import Alignment


##### ARGUMENT LIST ######

parser = argparse.ArgumentParser(description="Concatenates DNA data matrices")

# Main execution
main_exec = parser.add_argument_group("Main execution")
main_exec.add_argument("-in",dest="infile",nargs="+",help="Provide the input file name. If multiple files are provided, plase separated the names with spaces")
main_exec.add_argument("-if",dest="input_format",default="guess",choices=["fasta","nexus","phylip","guess"],help="Format of the input file(s). The default is 'guess' in which the program tries to guess the input format and genetic code automatically")
main_exec.add_argument("-of",dest="output_format",nargs="+",default="nexus",choices=["nexus","phylip","fasta","mcmctree"],help="Format of the ouput file(s). You may select multiple output formats simultaneously (default is '%(default)s')")
main_exec.add_argument("-o",dest="outfile",help="Name of the output file")

# Alternative modes
alternative = parser.add_argument_group("Alternative execution modes")
alternative.add_argument("-c",dest="conversion",action="store_const",const=True,help="Used for convertion of the input files passed as arguments with the -in option. This flag precludes the usage of the -o option, as the output file name is automatically generated based on the input file name")
#alternative.add_argument("-r",dest="reverse",help="Reverse a concatenated file into its original single locus alignments. A partition file similar to the one read by RAxML must be provided")
#alternative.add_argument("-z",dest="zorro",action="store_const",const=True,help="Use this option if you wish to concatenate auxiliary Zorro files associated with each alignment. Note that the auxiliary files must have the same prefix of the alignment file, with the addition of '_zorro.out'")
#alternative.add_argument("-zfile",dest="zorro_infile",nargs="*",default="_zorro.out",help="Provide the sufix for the concatenated zorro file (default is '%(default)s')")
#alternative.add_argument("-charset",dest="charset",help="Format the partition file of RAxML into a charset analogous to the Nexus block. A partition file similar to the one read by RAxML must be provided")
#alternative.add_argument("-partfile",dest="partfile",help="Format the charset analogous to the Nexus block into the partition file of RAxML")

# Formatting options
formatting = parser.add_argument_group("Formatting options")
formatting.add_argument("-model",dest="model_phy",default="WAG",choices=["DAYHOFF","DCMUT","JTT","MTREV","WAG","RTREV","CPREV","VT","BLOSUM62","MTMAM","LG"],help="This option only applies for the concatenation of protein data into phylip format. Specify the model for all partitions defined in the partition file")
formatting.add_argument("-interleave",dest="interleave",action="store_const",const="interleave",help="Specificy this option to write output files in interleave format (currently only supported for nexus files")
formatting.add_argument("-g",dest="gap",default="-",help="Symbol for gap (default is '%(default)s')")
formatting.add_argument("-m",dest="missing",default="n",help="Symbol for missing data (default is '%(default)s')")

# Data manipulation
manipulation = parser.add_argument_group("Data manipultation")
#manipulation.add_argument("-rm",dest="remove",nargs="*",help="Removes the specified taxa from the final alignment. Multiple taxa may be specified and separated by whitespace")
#manipulation.add_argument("--pickle-taxa",dest="pickle", choices=["dump","load"],help="Dump option: Only output a pickle object with the taxa names of the input alignment; Load option: loads the taxa names from a pickle object to be incorporated in the output alignment")

# Additional data checks
checks = parser.add_argument_group("Additional data checks")
checks.add_argument("--check",dest="check",action="store_const",const=True,help="Provides a final check for the lengths of the alignment sequences")

arg = parser.parse_args()

##### MAIN FUNCTIONS ######

def main_parser(alignment_list):
	""" Function with the main operations of ElConcatenero """
	
	# Defining main variables
	gap = arg.gap
	missing_sym = arg.missing
	conversion = arg.conversion
	input_format = arg.input_format
	output_format = arg.output_format
	outfile = arg.outfile
	interleave = arg.interleave
	model_phy = arg.model_phy

		# Defining output file name
	if arg.conversion == None and arg.outfile != None:
		outfile = "".join(arg.outfile)
	elif arg.conversion != None and arg.outfile != None:
		outfile = "".join(arg.outfile)
	elif arg.conversion != None and arg.outfile == None:
		if input_format in output_format:
			outfile = "".join(alignment_list).split(".")[0]+"_conv"
		else:
			outfile = "".join(alignment_list).split(".")[0]


	# Creating main instance of the parser
	# main_instance = ep.SeqUtils ()

	# The input file at this stage is not necessary
	# If just converting the partition file format do this and exit
	# if arg.charset != None:
	# 	partitions = main_instance.get_partitions(arg.charset)
	# 	output_instance = ep.writer()
	# 	output_instance.define_args(output_file = outfile)
	# 	output_instance.write_charset(partitions)
	# 	print ("Partitions formatting done!")
	# 	return 0
		
	# if arg.partfile != None:
	# 	partitions = main_instance.get_partitions(arg.partfile)
	# 	output_instance = ep.writer()
	# 	output_instance.define_args(output_file = outfile)
	# 	output_instance.write_partitions(partitions,model=arg.model_phy)
	# 	print ("Partitions formatting done!")
	# 	return 0

	# From here, the input file is mandatory
	# Guessing input format

	if len(alignment_list) == 1:

		alignment = Alignment.Alignment("".join(alignment_list))
		alignment.write_to_file(output_format, outfile)

	else:

		alignments = Alignment.AlignmentList(alignment_list)


	# if arg.input_format == "guess":
	# 	auto_format = ep.autofinder (alignment_list[0]) # The script will guess the input format of the first alignment in the list
	# 	if auto_format[0] != "unknown":
	# 		input_format = auto_format[0]
	# 		coding = auto_format[1][0]
	# 		missing_sym = auto_format[1][1]
	# 		print ("Input format set to %s.\nGenetic code set to %s." % (auto_format[0],auto_format[1][0]))
	# 	else:
	# 		print ("Input format could not be determined. Falling back to default values...")
	
	# # Parsing input file(s)
	# if len(alignment_list) == 1:
	# 	input_alignment = "".join(alignment_list)
	# 	alignment_storage = main_instance.read_alignment(input_alignment, input_format)
	# 	model = alignment_storage[3]
	# else:
	# 	alignment_storage = main_instance.read_alignments(alignment_list, input_format,missing=missing_sym)
	# 	model = alignment_storage[4]
	# 	if arg.zorro != None:
	# 		zorro_code = arg.zorro_infile
	# 		zorro_weigths = main_instance.zorro2rax(alignment_list,zorro_code)
	
	# If runnning reverse concatenation mode do this and quit
	# if arg.reverse != None:
	# 	partitions = main_instance.get_partitions(arg.reverse)
	# 	alignment_list = main_instance.reverse_concatenation(alignment_storage[0], partitions)
	# 	output_instance = ep.writer() # Initiating main output instance
	# 	output_instance.define_args(outfile, alignment_storage[1], coding, alignment_storage[2], alignment_storage[3],missing=missing_sym, models=model) # Providing main arguments for the output instance
	# 	output_instance.reverse_wrapper(alignment_list, "".join(output_format))
	# 	print ("Reverse concatenation complete!")
	# 	return 0

	# Removes specified taxa, if the option was declared. Otherwise, continue with the original alignment
	# if arg.remove != None:
	# 	alignment_dic, taxa_order = main_instance.rm_taxa(alignment_storage[0],arg.remove)
	# else:
	# 	alignment_dic = alignment_storage[0]
	# 	taxa_order = alignment_storage[1]
	
	# If the --pickle-taxa option is used, this code executes both loading and dumping operations
	# They can only be used separately, though.
	# if arg.pickle != None:
	# 	main_instance.pickle_taxa(alignment_dic,"".join(arg.pickle))
	# 	alignment_dic, taxa_order = main_instance.import_taxa(alignment_dic,len(list(alignment_dic.values())[0]),missing_sym)
		
	# Final consistency check of sequence lengths
	# if arg.check != None:
	# 	main_instance.check_sizes( alignment_dic, "Concatenated")					
	
					
	# Creating main output instance
	# output_instance = ep.writer()
	# output_instance.define_args(outfile, taxa_order, coding, alignment_storage[2], alignment_storage[3],missing=missing_sym, models=model)
	
	# Creating output file(s)
	# if "fasta" in output_format:
	# 	output_instance.fasta(alignment_dic,conversion)
	# if "phylip" in output_format:
	# 	output_instance.phylip(alignment_dic,model_phy,conversion)
	# if "nexus" in output_format:
	# 	output_instance.nexus(alignment_dic,conversion,interleave)
	# if arg.zorro != None:
	# 	output_instance.zorro(zorro_weigths)
	# if "mcmctree" in output_format:
	# 	if arg.charset == None and arg.partfile == None:
	# 		output_instance.mcmctree(alignment_dic, arg.partfile)
	# 	else:
	# 		if arg.partfile == None:
	# 			output_instance.mcmctree(alignment_dic, arg.charset)
	# 		else:
	# 			output_instance.mcmctree(alignment_dic, arg.partfile)

		
def main_check ():
	if arg.charset != None and arg.outfile == None:
		print ("ArgmentError: An output file must be provided with option '-o'\nExiting...")
		raise SystemExit
		
	if arg.charset != None or arg.partfile != None:
		return 0
		
	if arg.conversion == None and arg.outfile == None and arg.pickle == None and arg.reverse == None:
		print ("ArgumentError: If you wish to concatenate provide the output file name using the '-o' option. If you wish to convert a file, specify it using the '-c' option\nExiting...")
		raise SystemExit
		
	if len(arg.infile) == 1 and arg.conversion == None and arg.pickle == None and arg.reverse == None:
		print ("ArgumentError: Cannot perform concatenation of a single file. Please provide additional files to concatenate, or specify the conversion '-c' option.\nExiting...")
		raise SystemExit
		
	if arg.zorro != None and len(arg.infile) == 1:
		print ("ArgumentError: The '-z' option cannot be invoked when only a single input file is provided. This option is reserved for concatenation of multiple alignment files\nExiting...")
		raise SystemExit

	else:
		return 0
				
def main():
	main_check()
	if arg.conversion != None and len(arg.infile) > 1:
		for infile in arg.infile:
			main_parser([infile])		
	else:
		main_parser(arg.infile)

##### EXECUTION ######

main()
