#!/usr/bin/python3

# ElConcatenero v3.0.0-1
# Author: Diogo N Silva
# Last update: 15/10/2012
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
parser.add_argument("-of",dest="OutputFormat",nargs="+",default="nexus",choices=["nexus","phylip","fasta"],help="Format of the ouput file (default is '%(default)s')")
parser.add_argument("-c",dest="Conversion",action="store_const",const=True,help="Used for convertion of the input files passed as arguments with the -in option. This flag precludes the usage of the -o option, as the output file name is automatically generated based on the input file name.")
parser.add_argument("-o",dest="outfile",help="Name of the output file")
parser.add_argument("-z",dest="zorro",action="store_const",const=True,help="Use this option if you wish to concatenate auxiliary Zorro files associated with each alignment. Note that the auxiliary files must have the same prefix of the alignment file, with the addition of '_zorro.out'")
parser.add_argument("-zfile",dest="zorro_infile",nargs="*",default="_zorro.out",help="Please provide the sufix for the concatenated zorro file (default is '%(default)s')")
parser.add_argument("-in",dest="infile",nargs="+",required=True,help="Provide the input file name. If multiple files are provided, plase separated the names with spaces")
parser.add_argument("-code",dest="coding",default="DNA",choices=["DNA","Protein"],help="The coding of your molecular sequences")

arg = parser.parse_args()

# If you wish to include a whitespace (tab) between loci, set the variable to "yes", ortherwise "no".
tab_delimited_loci = "no"

##### MAIN FUNCTIONS ######

def main_parser(alignment_list):
	""" Function with the main operations of ElConcatenero """
	
	# Defining main variables
	gap = arg.gap
	missing = arg.missing
	coding = arg.coding
	conversion = arg.Conversion
	input_format = arg.InputFormat
	output_format = arg.OutputFormat
	output_file = arg.outfile

		
	# Creating main instance of the parser
	main_instance = ep.SeqUtils ()
	
	# Parsing input file(s)
	if len(alignment_list) == 1:
		input_alignment = "".join(alignment_list)
		alignment_storage = main_instance.read_alignment(input_alignment, input_format)
	else:
		alignment_storage = main_instance.read_alignments(alignment_list, input_format)
		if arg.zorro != None:
			zorro_code = arg.zorro_infile
			zorro_weigths = main_instance.zorro2rax(alignment_list,zorro_code)
	
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
	output_instance = main_instance.writer(output_file, alignment_storage[1], coding, alignment_storage[2], alignment_storage[3])
	
	# Creating output file(s)
	if "fasta" in output_format:
		output_instance.fasta(alignment_storage[0],conversion)
	if "phylip" in output_format:
		output_instance.phylip(alignment_storage[0],conversion)
	if "nexus" in output_format:
		output_instance.nexus(alignment_storage[0],conversion)
	if arg.zorro != None:
		output_instance.zorro(zorro_weigths)
		
def main_check ():
	if arg.Conversion == None and arg.outfile == None:
		print ("ArgumentError: If you wish to concatenate provide the output file name using the '-o' option. If you wish to convert a file, specify it using the '-c' option\nExiting...")
		raise SystemExit
	#~ if arg.Conversion != None and len(arg.infile) > 1:
		#~ print ("ArgumentError: Cannot convert multiple files implicitly. To convert multiple files independently use bash\nExample: for i in `ls`; do ElConcatenero [options]; done\nExiting...")
		#~ raise SystemExit
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
