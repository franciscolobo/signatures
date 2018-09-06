                        README for SIGNATURES
                     (last updated 23/08/2018)

AUTHORS
-=-=-=-

Francisco Pereira Lobo (franciscolobo@ufmg.br, franciscolobo@gmail.com)



1 - DESCRIPTION
-=-=-=-=-=-=-=-

SIGNATURES is a software that aims to compute and visualize genomic signatures
(usually a observed values of a DNA word - e.g., dinucleotides -  divided by 
the expected values of the same word for that sequence. SIGNATURES uses PERL
and R languages to extract genomic sequences from GenBank files, compute the
signatures and generate plots. For this reason users must install these
software to use SIGNATURES. Detailed installation instructions area available
in "INSTALL" file.



2 - HOW TO USE - IN BRIEF
-=-=-=-=-=-=-=-=-=-=-=-=-

There are *** steps to use SIGNATURES after downloading the source code:

A) Install SIGNATURES. For this purpose first install all POTION's dependencies
(third-party programs and perl modules) and set the parameters in 
'potion_config' file ('config_files' directory) according to the instructions
in "INSTALL" file.

(...)

C) 


3 - HOW TO USE - IN DETAIL
-=-=-=-=-=-=-=-=-=-=-=-=-=


3.1 - INSTALL
-=-=-=-=-=-=-

Install SIGNATURES following the instructions contained in 'INSTALL' file.
Check if everything went all right by executing SIGNATURES using the test
data distributed with it. Also, please ckeck the step-by-step guide available
at http://***
a detailed explanation of how to use SIGNATURES.

The install procedures will create a main directory called SIGNATURES-<version>.
This directory contains the following sub directories:

  bin          - contains a program (potion.pl) and module
                 (module_latest.pl), as well as a sub directory ('utils') that
                 contain several scripts to assist the initial steps needed to
                 obtain and format the sequence data to be analyzed with POTION.

  config_files - contains configuration files used to inform POTION about
                 the location of third-party software, as well as dummy empty
                 files used by POTION to create configuration files for third-
                 party progams and for POTION itself. Of special interest is
                 the file "potion_config", which users must edit in order to
                 tell POTION where it should look for third-party software.
  
  examples      - contains the two example datasets (TRYP and MYC) used to 
                 validate POTION. Use them to validate your local install
                 of POTION. For each dataset we provide the homology
                 relationship data, sequence data and POTION configuration
                 files to reproduce the results described in the article. 


3.2 - RUNNING SIGNATURES WITH EXAMPLE DATA
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

We distribute SIGNATURES with a example dataset ("***", in 'examples'
directory) that can be used to validate your install. To execute SIGNATURES
using these datasets edit the configuration file and point an output directory.
Please notice that SIGNATURES may produces a considerable amount of
intermediate files. A complete analysis of the *** dataset, for instance, can
easily create results of more than 5 Gb.


3.3 - SIGNATURE INPUT FILES
-=-=-=-=-=-=-=-=-=-=-=-=-=-

FORMATS
-=-=-=-

After installing SIGNATURES, users can start anayzing their own data. Each
"sample" is ideally a genbank file containing one or more sequences from a
genome. Users may also provide multi-fasta  files containing CDS data (for
codon-related analyses) or nucleotidic sequences of any nature (for the
remaining analyses).

Users must also provide a tab-delimited text file describing where sample files
are located and other sample-related metadata (e.g. genome strandedness, genome
group for clustering analyses).



SIGNATURE automatically parses each sample and reports several
nucleotide, dinucleotide and codon usage bias metrics for them in a visual way,
displaying several kinds of information.
reporting 


3.4 - ANALYZING YOR OWN DATA
-=-=-=-=-=-=-=-=-=-=-=-=-=-=

SIGNATURES is distributed with several scripts to automatize data acquisition and
allow users to analyze their own data with POTION ('bin/utils/' directory).
These scripts allow users to automatically download genbank text files and parse
them to extract CDS and protein data. To obtain homology data users can: 1) Use
the extracted protein data in OrthoMCL 1.4 and use its main results file or 2)
use a distributed script that parses OrthoXML data to the OrthoMCL results 
format in order to use data from specialized databases of homologous genes that
export homology information in OrthoXML format. A guide to use these scripts can
be found in the document "README_utils.txt", also located in 'utils' directory.


The behavior of POTION is dictated through a main parameter configuration file.
To create an empty configuration file to be filled with your experimental data
just execute POTION with the parameter '--create_conf'.


3.5 - SETTING UP SIGNATURE PARAMETERS
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

POTION contains two configuration files. One of them (install configuration
file) is used to configure POTION parameters that are not expected to change
after POTION install, such as the location of third-party programs, POTION's
root install directory, and the name of output files. The other configuration
file (project configuration file) contains parameters that are likely to change
each time POTION is executed, such as location of sequence files, filters to be
applied to data, or cutoffs to remove sequence and groups. Below follows the
explanation of all variables contained in these two files. Variables that point
to directory paths should end with "/" whenever possible. Comment lines in these
files are ignorated. A special variable ($potion_dir) is configured in install
configuration file, and should point to the complete path to directory
'POTION-<version>'. The variable potion_dir could be used within both POTION's
configuration files to refer to paths to other files/directories.
As an example, if you have:

potion_dir = /home/myname/POTION-<version>/

and want to store your results in a directory named "analysis_yyyymmdd", then
both of the following point to the same directory:

project_dir = /home/myname/POTION-<version>/projects/analysis_yyyymmdd
project_dir = $potion_dir/projects/analysis_yyyymmdd


3.5.1 - POTION install configuration file
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

This file is located in directory 'config_files/' and is called
'potion_config'. 

--------------
+ potion_dir +
--------------

potion_dir - path to the root directory of POTION.


---------------------------------
+ paths to third-party programs +
---------------------------------

codeml - path to the Codeml executable of PAML package
consense - path to the Consense executable of Phylip package
dnaml -  path to the Dnaml executable
muscle -  path to the MUSCLE executable
phipack - path to the Phi executable of PhiPack (http://www.maths.otago.ac.nz/~dbryant/software.html)
prank - path to the Prank executable
proml - path to the Proml executable of Phylip package
seqboot - path to the Seqboot executable of Phylip package
trimal - path to the Trimal executable
mafft - path to mafft executable
phyml - path to phyml executable

-------------------------
+ names of output files +
-------------------------

result_table - name of main table file. Default: final_results.out.

result_uncertain - name of fasta files containing sequences with some evidence
of positive selection, but not strong enough after FDR. Default: uncertain.out.

result_positive - name of files containing evidence of positive selection after
FDR. Default: positive.out

result_recombinants - name of files containing sequences with evidence of
recombination. Default: recombinants.out

summary - name of files containing summary of POTION execution. Default:
summary.out

tries - number of tries for each program, in case it crashes or returns an
error.


3.5.2 - POTION project configuration file
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

POTION contains a rich set of parameters to dictate how it should behave during
execution. Some of them are project-related, such as to define input and output
files and directories, or the number of processors to be used in the 
parallelized steps of POTION. Other parameters are the several filters and
customizations to taylor POTION to user-specific needs, such as to analyze only
putative 1-1 orthologs or paralogs, to remove sequence outliers in terms of
quality, length or identity. Finally, at several steps POTION allow users to
choose between distinct third-party software, or to configure some steps
of these software. Below follows the explanations of the parameters currently
available in POTION. An "M" indicates a mandatory parameter.

project parameters
-=-=-=-=-=-=-=-=-=

mode - POTION's main analysis mode. Currently POTION only supports "site" as
main analysis mode, future versions are expected to support "branch" and
"branch-site" modes as well.

CDS_dir_path - path to the directory containning the fasta CDS files for
each genome to be analyzed. (M)

homology_file_path - path to the file describing the homology relationships of
the CDS contained in fasta files, OrthoMCL 1.4 format. (M)

project_dir_path - path to the main directory where results will be created.
This directory will be created if it does not exist. (M)

max_processes - Defines the maximum number of processors to be used in the
parallelized steps of POTION algorithm. (M)

remove_identical - Defines if POTION should remove groups 100% identical
at nucleotide level. Possible values are "yes" or "no". (M)

verbose - 1 to print nice log messages telling you what is going on, 
0 otherwise. (M)

sequence/group parameters
-=-=-=-=-=-=-=-=-=-=-=-=-

groups_to_process - Defines which lines of the clusters file will be processed.
Use "all" to analyze every group, "-" to set groups between two given lines
(including the said lines). Use "!" to not process a specific line, can be used
with "-" to specify a set to not be processed. Useful if groups are taking too
long to finish. Use "," or ";" to set distinct line sets. Examples:
"1;4-10;12" will process groups described in lines 1, 4 to 10 and 12, "all;!3"
will process all groups except group in line 3, "all;!3-5" will process all
groups, except groups in lines 3 to 5. If not defined, POTION assumes "all".

behavior_about_bad_clusters - Defines what POTION should do if it finds a group
of homologous that contain any gene removed due to any filter. Possible options
are 0 (does not filter any sequence, not recommended), 1 (remove only sequences
and keep group) or 2 (remove all groups containning flagged sequences).
Sequences can be flagged by any of the following parameters below:
validation_criteria, absolute_min_sequence_size, absolute_max_sequence_size,
relative_min_sequence_size, relative_max_sequence_size, 
identity_most_similar_sequence, mean_sequence_identity. (M)

behavior_about_paralogs - Defines how POTION should handle mixed groups of
homologous genes that contain orthologs and/or paralogs. Possible values are
0 (analyze all sequences within group), 1 (remove all paralogs and analyze
only the remainning putative 1-1 orthologs), 2 (remove groups with paralogous),
3 (remove 1-1 orthologs and analyze all remainning paralogs together) and 4
(remove 1-1 orthologs and split remainning paralogs to analyze each species
individually). (M)

validation_criteria - Defines which quality criteria POTION should use to
remove individual sequences. Possible values are 0 (do not perform any sequence
quality filter), 1 (checks for valid start codons), 2 (checks for valid stop
codons), 3 (checks for sequence length multiple of three), 4 (checks for non
standard nucleotides) and "all" (applies all previous criteria).

additional_start_codons - these codons, plus the ones specified in codon table,
will be the valid start codons for validation purposes. (CTG, GTG), for
instance, will define these two codons as additional start codons. Parenthesis
and commas to enumerate more than one codon are mandatory.

additional_stop_codons: Same as above for start codons.

codon_table - Indicates the codon table id to be used by POTION for translation
and validation purposes, possible tables are defined in NCBI and bioperl.
Possible values are number from 1 to 22. See 
http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi for the definition
and taxonomic range of the 22 genetic codes currently supported. (M)

absolute_min_sequence_size - minimum absolute sequence length cutoff. Sequences
whose lengths are smaller than this value are removed from further analyses.
Possible values are integers between 0 and infinite.

absolute_max_sequence_size - maximum absolute sequence length cutoff. Sequences
whose lengths are greater than this value are removed from further analyses.
Possible values are integers between 0 and infinite.

relative_min_sequence_size - minimum relative sequence length cutoff. Sequences
smaller than mean|meadian times this value will be removed. Possible values are
numbers within the interval of 0 and 1. Mean/median values are calculated for
the whole group in order to allow the comparison.

relative_max_sequence_size - maximum relative sequence length cutoff. Sequences
greater than mean|meadian times this value will be removed. Possible values are
numbers between 1 and infinite. Mean/median values are calculated for the whole
group in order to allow the comparison.

sequence_size_average_metric - Dictates which metric will be calculated to
determine the minimum/maximum relative lengths ranges for  sequence removal
(configured in [relative_min|relative_max]_sequence_size). Possible values are
"mean" or "median". (M if relative_min or relative_max are defined).

min_group_identity - mean minimum group identity cutoff in pairwise sequence
alignments

max_group_identity - mean maximum group identity cutoff in pairwise sequence
alignments

group_identity_comparison - indicates the kind of sequence that will be used
when computing mean group identity. Possible values are "nt" or "aa"

min_sequence_identity - minimum (mean/median) sequence identity cutoff in
pairwise sequence alignments

max_sequence_identity - maximum (mean/median) sequence identity cutoff in
pairwise sequence alignemnts

sequence_identity_average_metric - dictates which metric will be calculated to
identify sequences with extreme values of identity (configured in 
[min|max]_sequence_identity). Possible values are "mean" and "median"

sequence_identity_comparison - indicates the kind of sequence that will be
used when computing sequence identity. Possible values are "nt" and "aa".

minimum_gene_number_per_cluster - lower bound cutoff for the minimum number of
genes in group after all filtering steps.

maximum_gene_number_per_cluster - upper bound cutoff for the minimum number of
genes in group after all filtering steps.

minimum_specie_number_per_cluster - lower bound cutoff for the minimum number 
of genes in group after all filtering steps.

maximum_specie_number_per_cluster - lower bound cutoff for the minimum number 
of genes in group after all filtering steps.

reference_genome_files - This is the final filter step applied by POTION. If
users chose a reference genome file POTION will only analyze groups that after
all previous filtering steps contain at least one gene from reference genome.
POTION will also return all the results found (sequence, IDs and coordinates)
in reference to the anchor genome. Possible value is one of the names of the 
fasta CDS files.


third-party software configuration
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

multiple_alignment - defines program used for multiple protein sequence 
alignment. Possible values are muscle|mafft|prank. Muscle is executed with
default parameters, prank is executed with the flags -twice (to run the 
analysis twice for each group) and -F (to correctly penalize the gaps) and 
mafft is executed with the flag --auto for auto-configuration. (M)

bootstrap - Number of bootstraps to be generated in phylogenetic analysis.
Possible values are integers. (M)

phylogenetic_tree - POTION currently supports proml|dnaml from phylip package 
and phyml. All of them calculate phyolgenetic trees using the maximum
likelihood approach on protein or nucleotide data, respectively. Possible 
values are proml|dnaml|phyml_nt|phyml_aa. (M)

phylogenetic_tree_speed - Both proml and dnaml contain a fast execution mode.
Use this parameter to turn it on or off. Possible values are fast|slow. (M)

rec_minimum_confirmations - This number dictates the minimum number of 
recombination tests with significant q-values (FDR-corrected p-values) in order
to infer recombination in a group. Possible values are 0, 1, 2 or 3. A value of
2, for instance, will demand significance in at least two recombination tests
to infer recombination occurrence.

rec_mandatory_tests - Dictates which recombination tests must be positive in
order to infer recombination. Possible values are phi, maxchi2 and nss. It is
possible to request more than one test. For instance, setting this variable as
phi maxchi2 will request POTION to find recombination for both phi and maxchi2
tests to classify a group as with evidence of recombination. The configuration
of this variable can be used with rec_minimum_confirmations. Setting these
variables for "phi" and "2", for instance, will request recombination
detection in at least two tests, and one of them must be phi. The definition of
the three recombination tests can be found in PhiPack package.

recombination_qvalue - Defines the q-value for recombination detection. It must
occur for all the specified tests. Possible values are numbers between 0 and 1.
A value of 0 will not look for recombination.

remove_gaps - This variable controls how POTION should trim columns of multiple
sequence alignments. If filled with numeric values between 0 and 1, POTION will
use this number as upper bound cutoff for the number of gaps in a given column
alignment. The values of "strict" or "strictplus" for this variable will
execute more complex column filters that take into account the conservation of
the neighborhood of a given column to trim sequence alignments. These filters 
are described in detail in the supplementary material of the article that 
describes trimal, the software that performs this step.

PAML_models - This variable configures which nested codon evolution models will
be computed by codeml. Possible values are m12, m78 and/or m8a8. To evaluate
all nested models write "m12 m78 m8a8" as value for this variable. (M)

pvalue - p-values for positive selection detection. Groups with p-values lower
than this cutoff and greater q-value cutoff are marked as "u" (uncertain) in
POTION's main result table and in the respective fasta files (please check the
respective sections below to obtain further information about POTION output
files).

qvalue - q-values for positive selection detection. Groups with q-values lower
than this cutoff are marked as "P" (positively-selected) in POTION's main 
results table and in the respective fasta files (please check the respective
sections below to obtain further information about POTION output files)


3.6 - POTION OUTPUT
-=-=-=-=-=-=-=-=-=-

POTION produces detailed intermediate and final results for groups evaluated.
All final results are stored in 'results' directory, intermediate files
are stored in 'intermediate_files' directory. Our software also produces 
detailed log information about each sequence and each group in, as well as
error log in files 'log' and 'log.err', respectively. POTION also create a copy
of the configuration file used to produce a given analysis.


main output files
-=-=-=-=-=-=-=-=-

Currently our software returns as final results the following files, whose
names can be configured in 'potion_config' file ('config_files' directory):


--------------------
+ final_results.out +
--------------------

Tabular flat file containing the statistical data and final classification for
each valid group (not removed due to any filter). The values contained in this
file are:

 - cluster name;

 - final classification for nested models M1a2 and M78 (possible values are "P"
   for groups with q-values smaller than cutoff, "u", for groups with q-values
   greater than cutoff and p-values smaller than cutoff, "n", for negative
   groups, and "-", fro groups where it was not possible to compute p-values);

 - -lnL values for models evaluated (M1a, M2, M7, M8a, M8). Possible values are
   negative number or "-" symbols for groups where it was not possible to 
   compute -lnL values;

 - D values, the difference of likelihoods for neutral models (M1a, M7 and M8a)
  when compared with models that also allow sites under positive selection (M2,
  M8 and M8, respectively).

 - p-values computed for each pair of nested models evaluated. Calculated
   through a Likelihood Ratio Test (LRT).

 - q-values computed for each pair of nested models evaluated. Calculated
   through a False Discovery Rate correction.


---------------------------
+ allvalid_for_enrichment +
---------------------------

Fasta flat file containing all valid groups after filtering steps, but prior to
recombination detection. To be used for enrichment analysis as background 
distribution if users want to compare positively selected genes against the
background of valid genes without considering recombination.


-----------------------------------
+ non_recombinants_for_enrichment +
-----------------------------------

Fasta flat file containing all valid groups after filtering steps after 
recombination detection. To be used for enrichment analysis as background
distribution if users want to compare positively selected genes against the
background of valid genes with no evidence of recombination.


----------------
+ positive.out +
----------------

Flat file that contains information about sites under positive selection found
in each group of homologous genes. Each group is described in four lines:

First line:
<Group ID>: <gene ID 1> <gene ID 2> ... <gene ID n>

Second line:
><gene id> (model where positive selection was detected)

Third line:
sequence

Fourth line:
sites under positive selection, where "-" indicates no evidence, "+" indicates
Bayes Empyrical Bayes (BEB) values greater than 0.95, "*" indicates BEB values
greater than 0.99.


----------------------------
+ positive.out.interleaved +
----------------------------

Same as positive.out, but sequence and sites under positive selection are 
displayed in interleaved format to make it easier to visualize where are
the sites under positive selection.


-----------------------------------------
+ positive.out.[aa|nt].*_for_enrichment +
-----------------------------------------

Fasta flat files containing protein or nucleotide sequence data for groups of
homologous found as positively selected (q-values smaller than cutoff). POTION
will create at most three distinct files to store genes found as positively
selected for m2, m8 and both pairs of nested models.


-----------------------------------------
+ recombinants.[aa|nt].*_for_enrichment +
-----------------------------------------

Fasta flat files containing protein or nucleotide sequence data for groups of
homologous with evidence of recombination (q-values smaller than cutoff).


---------------
+ summary.out +
---------------

Flat file containing information about POTION execution, such as total number
of valid clusters, time spent, theoretical total time spent if POTION was
executed with a single processor, total time spent in each of the 
parallelized steps of POTION.


---------------
+ uncertain.* +
---------------

Flat files containing similar information as the one produced for positively
selected genes. Naming conventions are the same.



intermediate files
-=-=-=-=-=-=-=-=-=

The directory 'intermediate_files' contain a directory for each valid group
submitted to the parallelized portion of POTION (starting with multiple protein
sequence alignment). Inside each group directory are all intermediate files
created during POTION execution. File contents are as follows:


id_names - links ids to internal ids used by POTION. Needed since several
programs/formats have requests that may not be fulfilled by initial gene ids
such as phylip maximum sequence name length of ten.

ORTHOMCL123.group_status - if exists, this file contains information about the
reasons why the computation of this group could not go trough the entire POTION
pipeline. It may be the error message of third-party software, if POTION was
able to capture it, or some group filtering criteria.

ORTHOMCL123.cluster.aa.fa - amino acid fasta sequence file 

ORTHOMCL123.cluster.nt.fa - nucleotide fasta sequence file

ORTHOMCL123.cluster.aa.fa.aln - aligned amino acid fasta sequence file

ORTHOMCL123.cluster.aa.fa.aln.1.fas - intermediate alignment file

ORTHOMCL123.cluster.aa.fa.aln.[aa|nt].phy - aligned sequence in phylip format
for phylogenetic analysis.

ORTHOMCL123.cluster.aa.fa.aln.[aa|nt].phy.trim - aligned trimmed sequence in
phylip format for phylogenetic analysis.

ORTHOMCL123.cluster.aa.fa.aln.[aa|nt].phy.trim.boot - bootstrap file in phylip
format for phylogenetic analysis

ORTHOMCL123.cluster.aa.fa.aln.[aa|nt].phy.trim.boot.log - execution log for
bootstrap program, phylip package

ORTHOMCL123.cluster.aa.fa.aln.aa.phy.trim.conf.boot - configuration file used
to create bootstrap file using bootstrap program, phylip package

ORTHOMCL123.cluster.aa.fa.aln.aa.phy.trim.conf.final_tree - configuration file
used to create consensus tree using consensus program, phylip package

ORTHOMCL123.cluster.aa.fa.aln.aa.phy.trim.conf.tree - configuration file used
to create phylogenetic trees using program [dna|pro]ml, phylip package, or
PhyML software

ORTHOMCL123.cluster.aa.fa.aln.aa.phy.trim.final_tree - consensus tree, newick
format, produced using consensus tree, phylip package

ORTHOMCL123.cluster.aa.fa.aln.aa.phy.trim.final_tree.log - execution log of
consensus program, phylip package

ORTHOMCL123.cluster.aa.fa.aln.aa.phy.trim.tree - trees produced using bootstrap
files, produced using [dna|pro]ml, phylip package or PhyML software

ORTHOMCL123.cluster.aa.fa.aln.aa.phy.trim.tree.log - execution log of
[dna|pro]ml, phylip package or PhyML software

ORTHOMCL123.cluster.aa.fa.aln.nt.fa - codon alignment data, produced using
protein alignemnt data

ORTHOMCL123.cluster.aa.fa.aln.nt.trim.fa - codon alignment trimmed data,
produced using trimal

ORTHOMCL123.cluster.aa.fa.aln.nt.phy - codon alignment data, phylip format

ORTHOMCL123.cluster.aa.fa.aln.nt.phy.trim - trimmed codon alignment data,
phylip format

ORTHOMCL123.cluster.aa.fa.aln.nt.phy.trim.paml.model[1|2|7|8] - result of
codeml execution for model M1a, M2, M7 or M8.

ORTHOMCL123.cluster.aa.fa.aln.nt.phy.trim.paml.model[1|2|7|8].config - 
configuration file used to execute codeml for model M1a, M2, M7 or M8.

ORTHOMCL123.cluster.aa.fa.aln.nt.phy.trim.paml.model[1|2|7|8].log - execution
log of codeml for model M1a, M2, M7 or M8.

ORTHOMCL123.cluster.aa.fa.aln.pairs - contains information about group and 
sequence identity metrics. Calculated using trimal.

ORTHOMCL123.trimal_cols.aa - contains columns positions filtered out by trimal.
Used to trim codon alignment data.

Phi, Phi.inf.list, Phi.inf.sites, Phi.poly.unambig.sites - result and log files
produced by phipack package when searching for recombination

2NG.dN, 2NG.dS, 2NG.t, 4fold.nuc, lnf, rst, rst1, rub - intermediate files
produced by codeml

time_f_tree, time_id_rec, time_model_1, time_model_2, time_model_7, 
time_model_8 - time used to compute multiple sequence alignment and phylogenetic
tree, recombination, and codeml models 1, 2, 7 and 8


log files
-=-=-=-=-

POTION generates an execution log where detailed information about each group
and sequence can be found, including the reasons for removal from further
analyses.


3.7 - FREQUENTLY ASKED QUESTIONS
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Q - How did you choose the parameters for the "soft", "standard" and "hard"
filters?

A - Intuition and a lot of manual inspection of alignment data! In brief,
"soft" may contain some false-positives, but is useful if you did not find any
result when using POTION in "standard" or "hard" mode and wish to see if you
find anything suspicious for further investigation. "standard" appears to give
good results and removes most of the noisy data, and usually is our first
choice when analysing data with POTION for the first time. "hard" is just what
it means to be: it will remove almost all data with looks barely suspicious,
but what remains is very likely to represent true positive selection, at least
in our current knowledge of POTION's behavior.

Q - How POTION behaves if there is more than one gene from the anchor genome
in a given group?

A - POTION will report only the longest ORF if you are running any analysis
where more than one gene from the same species is to be reported. Possible
cases are analyses of paralogs or analyses anchored to a given genome.


3.8 - DEVELOPMENT NOTES/KNOWN ISSUES
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Some notes about the third-party programs/module:


- New versions of BioPerl could cause POTION to get stuck when transforming
 alignment files from one format to another. Newer versions of POTION (from 
 1.1.2 on) do not rely on AlignIO BioPerl module for this step.

- Prank currently uses the parameters -twice, -F and -quiet; the second file
  produced as output is used as the alignment

- Some programs can crash during execution even if they would run normally if 
  done manually. Programs that have any bug regarding memory leaks have a
  higher likelihood to run into it when used in a highly-parallelized pipeline;
  POTION will attempt to isolate the problem and re-run the program at least
  three times (can be set by the user) before excluding the group from the
  analysis.

- Prank may abort a few times for a large number of processes (we are not sure,
  but it appears to be due to memory issues; Potion will attempt to
  reinitialize it until it runs properly. If rerunning doesn't solve it, we
  recommend running with less processes, reducing the parameter 'max_processes'
  until the problem stops occurring.

- If you want to halt prank during a run, use 'killall perl' followed by
  'killall prank'. However, it doesn't prevent possible errors due to memory
  overwrite.

- "version `GLIBC_x.xx' not found" is a possible message when running PAML,
  Phylip executables and other programs written in C. This indicates that you
  compiled them with a newer C library than the one in the computer you are
  using. Update your C library (GLIBC) from
  http://www.gnu.org/software/libc/index.html or recompile the executable with
  your current library.
