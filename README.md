## About PhD_Easy

Note: Name under consideration

PhD_Easy is a python program developed for the UNIX OS (but should work in any operating system with python 3.x) that converts, concatenates and manipulates common population genetics and phylogenetics data file types (such as Fasta, Phylip and Nexus). For the moment, the program works through command line exclusively, which makes it fast and scriptable, but a GUI version is being worked on. It requires the wingman library.

It currently supports as input file formats:

- Fasta
- Phylip
- Nexus (Alignment files and partition files)
- RAxML partition files

And its able to convert/concatenate any of these file formats into:

- Fasta
- Phylip (Including the *.parts file required for multilocus datasets in RAxML)
- Nexus
- ZORRO weigths files (these weigths files that are associated with their respective alignment files can be jointly concatenated)

Recently, PhD_Easy has received several new features:

- Collapse alignments into unique haplotypes. Taxa in the output alignment will be Haplotypes and their correspondance to the original taxa will be provided in a secondary file with the '.haplotypes' extension
- Code gaps into a binary state matrix that is appended at the end of the alignment, following the method of Simmons and Ochoterena (2000). This will only work when the output format is Nexus, and it can be red by MrBayes
- Columns of the alignment with excessive missing data can be filtered according to two user-specified thresholds (a threshold for the maximum proportion of gap and missing data, and another for missing data only).

Other minor improvements include:

- Specific taxa can be removed from the output alignments when converting/concatenating files
- In the Nexus output format, the outgroup may be set when converting/cocatenating 

There is no need for instalation. The only requirement is that the wingman library must be on the same directory as PhD_Easy.py (or you can use any other way to let the main script know where the module is). I do recommend, to make it easier to call the program, that you add it to your $PATH variable. This can be done by declaring "export $PATH=$PATH:/path/to/PhD\_Easy" on your .bashrc, or your shell specific rc file.

Finally, please note that PhD\_Easy.py is not immune to bugs, and I'll be happy to know about them through my e-mail (o.diogosilva@gmail.com) so that I can fix them. Any suggestions or comments are also welcome.

### Options

PhD_Easy.py has the following options (which can also be consulted by typing "PhD_Easy.py -h" or "PhD_Easy.py --help" in the command line):

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

  -r *REVERSE*            **Reverse a concatenated file into its original single
                        locus alignments. A partition file similar to the one
                        read by RAxML or a file containing a NEXUS charset block must be provided**

  -z ZORRO, --zorro-suffix *ZORRO*
                        **Use this option if you wish to concatenate auxiliary
                        Zorro files associated with each alignment. Provide
                        the sufix for the concatenated zorro file**

  -p PARTITION_FILE, --partition-file *PARTITION_FILE*
                        **Using this option and providing the partition file
                        will convert it between a RAxML or Nexus format**

  -collapse            **Use this flag if you would like to collapse the input
                        alignment(s) into unique haplotypes**

  -gcoder               **Use this flag to code the gaps of the alignment into a
                        binary state matrix that is appended to the end of the
                        alignment**

  -filter *FILTER FILTER*
                        **Use this option if you wish to filter the alignment's
                        missing data. Along with this option provide the
                        threshold percentages for gap and missing data,
                        respectively (e.g. -filter 50 75 - filters alignments
                        columns with more than 50% of gap+missing data and
                        columns with more than 75% of true missing data)**



####Formatting options:

   -model *{DAYHOFF,DCMUT,JTT,MTREV,WAG,RTREV,CPREV,VT,BLOSUM62,MTMAM,LG}*
                        **This option only applies for the concatenation of
                        protein data into phylip format. Specify the model for
                        all partitions defined in the partition file (default
                        is 'LG')**

  -interleave           **Specificy this option to write output files in
                        interleave format (currently only supported for nexus
                        files**

  -g *GAP*                **Symbol for gap (default is '-')**

  -m *MISSING*            **Symbol for missing data (default is 'n')**


####Data manipultation:

  -rm *[REMOVE [REMOVE ...]]*
                        **Removes the specified taxa from the final alignment.
                        Multiple taxa may be specified and separated by
                        whitespace**
                        
  -outgroup *[OUTGROUP_TAXA [OUTGROUP_TAXA ...]]*
                        **Provide taxon names/number for the outgroup (This
                        option is only supported for NEXUS output format
                        files)**xxxxxxxx

####Miscellaneous:

  -quiet                Removes all terminal output

#####Note: The order of the options does not matter.
		
### Usage examples

##### Conversion (Fasta to Nexus)

PhD_Easy.py -c -of nexus -in input_file.fas

##### Conversion (Phylip to Fasta)

PhD_Easy.py -c -of nexus -in input_file.phy

##### Conversion (Nexus to Phylip and fasta)

PhD_Easy.py -c -of phylip fasta -in input_file.nex

##### Conversion (Multiple Fasta files to Nexus)

PhD_Easy.py -c -of nexus -in input_file1.fas input_file2.fas input_file3.fas input_file4.fas [...] input_fileN.fas

or

PhD_Easy.py -c -of nexus -in *.fas

##### Concatenation (Nexus to Nexus)

ElConcetenero.py -of nexus -in input_file1.nex input_file2.nex input_file3.nex (...) input_fileN.nex -o concatenated_file

or

ElConcetenero.py -of nexus -in *.nex -o concatenated_file

##### Concatenation (Fasta to Phylip)

PhD_Easy.py -of phylip -in input_file1.fas input_file2.fas (...) input_fileN.fas -o concatenated_file

or

PhD_Easy.py -of phylip -in *.fas -o concatenated_file

#####Collapse alignment 

PhD_Easy.py -in *.fas -of fasta -collapse -o Collapsed_alignment

##### Coding gaps

PhD_Easy.py -in *.fas -of nexus -gcoder -c 

##### Filtering missing data

PhD_Easy.py -in *.fas -of fasta nexus phylip -filter 50 75

##### Remove taxa

PhD_Easy.py -in *.fas -of fasta -rm taxon1 taxon2 taxon3 (...) taxonN

#### ToDo

- Add outgroup specification for nexus format