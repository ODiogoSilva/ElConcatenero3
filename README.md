## About ElConcatenero3

ElConcatenero3 is a python program developed for the UNIX OS (but should work in any operating system with python 3.x) that converts and concatenates common population genetics and phylogenetics data file types (such as Fasta, Phylip and Nexus). The program works through command line exclusively, which makes it fast and scriptable, a usefull feature when dealing with lots of files simultaneously. It uses the ElParsito3.py module, which is where the file parsing and writting classes are.

It currently supports as input file formats:

- Fasta
- Phylip
- Nexus

And its able to convert/concatenate any of these file formats into:

- Fasta
- Phylip (Including the *.parts file required for multilocus datasets in RAxML)
- Nexus
- ZORRO weigths files (these weigths files that are associated with their respective alignment files can be jointly concatenated)

There is no need for instalation. The only requirement is that the ElParsito3.py module must be on the same directory as ElConcatenero3.py (or you can use any other way to let the main script know where the module is). I do recommend, to make it easier to call the program, that you add it to your $PATH variable. This can be done by declaring "export $PATH=$PATH:/path/to/ElConcatenero3" on your .bashrc, or your shell specific rc file.

Finally, please note that ElConcatero3.py is not immune to bugs, and I'll be happy to know about them through my e-mail (o.diogosilva@gmail.com) so that I can fix them. Any suggestions or comments are also welcome.

### Options

ElConcatenero3.py has the following options (which can also be consulted by typing "ElConcatenero3.py -h" or "ElConcatenero3.py --help in the command line):

####optional arguments:

  -h, --help            show this help message and exit

####Main execution:

  -in *INFILE [INFILE ...]*
						**Provide the input file name. If multiple files are
                        provided, plase separated the names with spaces**
                        
  -if *{fasta,nexus,phylip,guess}*
                        **Format of the input file(s). The default is 'guess' in
                        which the program tries to guess the input format and
                        genetic code automatically**
                        
  -of *{nexus,phylip,fasta} [{nexus,phylip,fasta} ...]*
                        **Format of the ouput file(s). You may select multiple
                        output formats simultaneously (default is 'nexus')**
                        
  -o *OUTFILE*           **Name of the output file**

####Alternative execution modes:

  -c                    **Used for convertion of the input files passed as
                        arguments with the -in option. This flag precludes the
                        usage of the -o option, as the output file name is
                        automatically generated based on the input file name**
                        
  -r REVERSE            **Reverse a concatenated file into its original single
                        locus alignments. A partition file similar to the one
                        read by RAxML must be provided**
                        
  -z                    **Use this option if you wish to concatenate auxiliary
                        Zorro files associated with each alignment. Note that
                        the auxiliary files must have the same prefix of the
                        alignment file, with the addition of '_zorro.out'**
                        
  -zfile *[ZORRO_INFILE [ZORRO_INFILE ...]]*
                        **Provide the sufix for the concatenated zorro file
                        (default is '_zorro.out')**
                        
  -charset *CHARSET*     **Format the partition file of RAxML into a charset
                        analogous to the Nexus block. A partition file similar
                        to the one read by RAxML must be provided**
                        
  -partfile *PARTFILE*    **Format the charset analogous to the Nexus block into
                        the partition file of RAxML**

####Formatting options:

  -model *{DAYHOFF,DCMUT,JTT,MTREV,WAG,RTREV,CPREV,VT,BLOSUM62,MTMAM,LG}*
                        **This option only applies for the concatenation of
                        protein data into phylip format. Specify the model for
                        all partitions defined in the partition file**
                        
  -interleave           **Specificy this option to write output files in
                        interleave format (currently only supported for nexus
                        files**
                        
  -g GAP                *Symbol for gap (default is '-')*
  
  -m MISSING            *Symbol for missing data (default is 'n')*

####Data manipultation:

  -rm *[REMOVE [REMOVE ...]]*
                        **Removes the specified taxa from the final alignment.
                        Multiple taxa may be specified and separated by
                        whitespace**
                        
  --pickle-taxa *{dump,load}*
                        **Dump option: Only output a pickle object with the taxa
                        names of the input alignment; Load option: loads the
                        taxa names from a pickle object to be incorporated in
                        the output alignment**

####Additional data checks:

  --check               **Provides a final check for the lengths of the
                        alignment sequences**

#####Note: The order of the options does not matter.
		
### Usage

##### Conversion (Fasta to Nexus)

ElConcatenero3.py -c -of nexus -in input_file.fas

##### Conversion (Phylip to Fasta)

ElConcatenero3.py -c -of nexus -in input_file.phy

##### Conversion (Nexus to Phylip and fasta)

ElConcatenero3.py -c -of phylip fasta -in input_file.nex

##### Conversion (Multiple Fasta files to Nexus)

ElConcatenero3.py -c -of nexus -in input_file1.fas input_file2.fas input_file3.fas input_file4.fas [...] input_fileN.fas

or

ElConcatenero3.py -c -of nexus -in *.fas

##### Concatenation (Nexus to Nexus)

ElConcetenero.py -of nexus -in input_file1.nex input_file2.nex input_file3.nex (...) input_fileN.nex -o concatenated_file

or

ElConcetenero.py -of nexus -in *.nex -o concatenated_file

##### Concatenation (Fasta to Phylip)

ElConcatenero3.py -of phylip -in input_file1.fas input_file2.fas (...) input_fileN.fas -o concatenated_file

or

ElConcatenero3.py -of phylip -in *.fas -o concatenated_file

#### ToDo

- Ability to parse concatenated datasets, either to separate them or to further concatenate
