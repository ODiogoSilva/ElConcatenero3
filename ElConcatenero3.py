#!/usr/bin/python3

# ElConcatenero v3.0.1
# Author: Diogo N Silva
# Last update: 25/10/2012
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
import ElParsito3 as ep


##### ARGUMENT LIST ######

parser = argparse.ArgumentParser(description="Concatenates DNA data matrices")
parser.add_argument("-g",dest="gap",default="-",help="Symbol for gap (default is '%(default)s')")
parser.add_argument("-m",dest="missing",default="n",help="Symbol for missing data (default is '%(default)s')")
parser.add_argument("-if",dest="InputFormat",default="fasta",choices=["fasta","nexus","phylip"],help="Format of the input file(s) (default is '%(default)s')")
parser.add_argument("-of",dest="OutputFormat",nargs="+",default="nexus",choices=["nexus","phylip","fasta"],help="Format of the ouput file(s). You may select multiple output formats simultaneously (default is '%(default)s')")
parser.add_argument("-c",dest="Conversion",action="store_const",const=True,help="Used for convertion of the input files passed as arguments with the -in option. This flag precludes the usage of the -o option, as the output file name is automatically generated based on the input file name.")
parser.add_argument("-o",dest="outfile",help="Name of the output file")
parser.add_argument("-z",dest="zorro",action="store_const",const=True,help="Use this option if you wish to concatenate auxiliary Zorro files associated with each alignment. Note that the auxiliary files must have the same prefix of the alignment file, with the addition of '_zorro.out'")
parser.add_argument("-zfile",dest="zorro_infile",nargs="*",default="_zorro.out",help="Provide the sufix for the concatenated zorro file (default is '%(default)s')")
parser.add_argument("-in",dest="infile",nargs="+",required=True,help="Provide the input file name. If multiple files are provided, plase separated the names with spaces")
parser.add_argument("-code",dest="coding",default="DNA",choices=["DNA","Protein"],help="The coding of your molecular sequences")
parser.add_argument("-rm",dest="remove",nargs="*",help="Removes the specified taxa from the final alignment. Multiple taxa mais be specified and separated by whitespace")
parser.add_argument("--pickle-taxa",dest="pickle", choices=["dump","load"],help="dump option: Only output a picke object with the taxa names of the input alignment; load option: loads the taxa names from a pickle object to be incorporated in the output alignment")
parser.add_argument("--check",dest="check",action="store_const",const=True,help="Provides a final check for the lengths of the alignment sequences")

arg = parser.parse_args()

##### MAIN FUNCTIONS ######

def main_parser(alignment_list):
	""" Function with the main operations of ElConcatenero """
	
	# Defining main variables
	gap = arg.gap
	missing_sym = arg.missing
	coding = arg.coding
	conversion = arg.Conversion
	input_format = arg.InputFormat
	output_format = arg.OutputFormat
	output_file = arg.outfile

		
	# Creating main instance of the parser
	main_instance = ep.SeqUtils (missing_sym)
	
	# Parsing input file(s)
	if len(alignment_list) == 1:
		input_alignment = "".join(alignment_list)
		alignment_storage = main_instance.read_alignment(input_alignment, input_format)
	else:
		alignment_storage = main_instance.read_alignments(alignment_list, input_format)
		if arg.zorro != None:
			zorro_code = arg.zorro_infile
			zorro_weigths = main_instance.zorro2rax(alignment_list,zorro_code)
	
	# Removes specified taxa, if the option was declared. Otherwise, continue with the original alignment
	if arg.remove != None:
		alignment_dic, taxa_order = main_instance.rm_taxa(alignment_storage[0],arg.remove)
	else:
		alignment_dic = alignment_storage[0]
		taxa_order = alignment_storage[1]
	
	# If the --pickle-taxa option is used, this code executes both loading and dumping operations
	# They can only be used separately, though.
	if arg.pickle != None:
		main_instance.pickle_taxa(alignment_dic,"".join(arg.pickle))
		alignment_dic, taxa_order = main_instance.import_taxa(alignment_dic)
		
	# Final consistency check of sequence lengths
	if arg.check != None:
		main_instance.check_sizes( alignment_dic, "Concatenated")
				
	# Defining output file name
	if arg.Conversion == None and arg.outfile != None:
		output_file = "".join(arg.outfile)
	elif arg.Conversion != None and arg.outfile != None:
		output_file = "".join(arg.outfile)
	elif arg.Conversion != None and arg.outfile == None:
		if input_format in output_format:
			output_file = "".join(alignment_list).split(".")[0]+"_conv"
		else:
			output_file = "".join(alignment_list).split(".")[0]
					
	# Creating main output instance
	output_instance = ep.writer(output_file, taxa_order, coding, alignment_storage[2], alignment_storage[3],missing=missing_sym)
	
	# Creating output file(s)
	if "fasta" in output_format:
		output_instance.fasta(alignment_dic,conversion)
	if "phylip" in output_format:
		output_instance.phylip(alignment_dic,conversion)
	if "nexus" in output_format:
		output_instance.nexus(alignment_dic,conversion)
	if arg.zorro != None:
		output_instance.zorro(zorro_weigths)
		
def main_check ():
	if arg.Conversion == None and arg.outfile == None and arg.pickle == None:
		print ("ArgumentError: If you wish to concatenate provide the output file name using the '-o' option. If you wish to convert a file, specify it using the '-c' option\nExiting...")
		raise SystemExit
	if len(arg.infile) == 1 and arg.Conversion == None and arg.pickle == None:
		print ("ArgumentError: Cannot perform concatenation of a single file. Please provide additional files to concatenate, or specify the conversion '-c' option.\nExiting...")
		raise SystemExit
	if arg.zorro != None and len(arg.infile) == 1:
		print ("ArgumentError: The '-z' option cannot be invoked when only a single input file is provided. This option is reserved for concatenation of multiple alignment files\nExiting...")
		raise SystemExit
				
def main():
	main_check()
	if arg.Conversion != None and len(arg.infile) > 1:
		for infile in arg.infile:
			main_parser([infile])		
	else:
		main_parser(arg.infile)

##### EXECUTION ######

main()
