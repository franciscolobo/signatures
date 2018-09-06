#!/usr/bin/perl

################################################################################
##                                                                            ##
## Copyright 2018 - Universidade Federal de Minas Gerais                      ##
## Authors: Tarcisio Jose Domingos Coutinho/Francisco Pereira Lobo            ##
## this program is free software: you can redistribute it and/or modify       ##
## it under the terms of the GNU General Public License as published by the   ##
## Free Software Foundation, version 3 of the License.                        ##
##                                                                            ##
## signatures.pl is distributed in the hope that it will be useful,           ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of             ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                       ##
## See the GNU General Public License for more details.                       ##
##                                                                            ##
## You should have received a copy of the GNU General Public License          ##
## along with rscu.pl (file: COPYING).                                        ##
##                                                                            ##
## If not, see <http://www.gnu.org/licenses/>.                                ##
##                                                                            ##
################################################################################

#This is the main script for SIGNATURES. It works as a wrapper to execute
#SIGNATURES

use strict;
use warnings;

my $version = "1.0.1";

# if users do not provide a configuration file or ask for help
if ((!$ARGV[0]) || ($ARGV[0] eq "-h") || ($ARGV[0] eq "--help")) {
  usage($version);
}

#path to configuration file
my $conf_file_path = $ARGV[0];

#stores parameters needed for execution
my %pars;

#stores genome-related data;
my %data; 

#stores paths to third-party software
my %path2external;

#parsing configuration file
print("Parsing configuration file ($conf_file_path)...\n");
my $tmp_hash = parse_config_file($conf_file_path);
%pars = %{$tmp_hash};
print("Done.\n\n");

#parsing metadata file
print("Parsing metadata file ($pars{metadata_file_path})...\n");
$tmp_hash = "";
$tmp_hash = parse_metadata_file($pars{metadata_file_path});
%data = %{$tmp_hash};
print("Done.\n\n");

#parsing paths to third-party software
print("Parsing path to external executables...\n");
$tmp_hash = "";
my $path2ext = "config_files/path2external.tab";
$tmp_hash = parse_path2external_file($path2ext);
%path2external = %{$tmp_hash};

print "SIGNATURES found the following third-party software:\n\n";
foreach my $key (keys %path2external) {
  print("$key\t$path2external{$key}\n");
}
#my $a = <STDIN>;

print("Validating execution parameters...\n");
#checking if parameters are what they should be
validate_parameters(\%pars);
print("Done.\n\n");

#checking if metadata is what it should be. Specifically:
## checking if file exists
## checking if all groups have at least two genomes
#validate_parameters()

#creating directory structure for output
print("Creating output directory structure...\n");
create_output_dir_structure($pars{outdir_path});
print("Done.\n\n");

#clean data
print("Parsing and cleaning sequence data (this may take a while...)\n");
my $task = "clean";
compute_task("clean", \%pars, \%data);
print("Done.\n\n");

#computing sequence individual tasks
my @tasks = split(/\s*,\s*/, $pars{tasks});
foreach my $task (@tasks) {
  print ("Computing task $task...\n");
  compute_task($task, \%pars, \%data);
  print("Done.\n\n");
}

#creating output report
print("Creating output report...\n");
create_output_report(\@tasks, \%pars, \%data);
print("Done.\n\n");

#subroutines
sub create_output_dir_structure {
  my $root_dir = $_[0];
  my $command = "";
  my $curr_dir = ""; #dir to be created;
  $curr_dir = $root_dir."/tmp/";
  if (-e $curr_dir) { #checking if directory exists, creating if not
  
  } else {
    $command = "mkdir ".$root_dir."/tmp/";
    system($command);
  }

  $curr_dir = $root_dir."/tmp/nuc_per_pos_files/";
  if (-e $curr_dir) { #checking if directory exists, creating if not
    
  } else {
    $command = "mkdir ".$root_dir."/tmp/parsed_files/";
    system($command);
    $command = "mkdir ".$root_dir."/tmp/valid_files/";
    system($command);
    $command = "mkdir ".$root_dir."/tmp/dinuc_files/";
    system($command);
    $command = "mkdir ".$root_dir."/tmp/RSCU_files/";
    system($command);
    $command = "mkdir ".$root_dir."/tmp/ENC_files/";
    system($command);
    $command = "mkdir ".$root_dir."/plots/";
    system($command);
    $command = "mkdir ".$root_dir."/logs/";
    system($command);
  $command = "mkdir ".$root_dir."/tmp/nuc_per_pos_files/";
  system($command);
  }
}

sub compute_task {
  my $task = $_[0];
  my %pars = %{$_[1]};
  my %data = %{$_[2]};

  my @supported_tasks = ("RSCU", "DINUC_ODDS_RATIO", "ENC", "NUC_PER_POSITION");

  if ($task =~ /^clean$/i) {
    print("\t* cleaning sequence data\n");
    clean_sequence_data(\%pars, \%data);
  }

  elsif ($task =~ /^RSCU$/i) {

    print("\t* computing RSCU values\n");
    compute_RSCU(\%pars, \%data);

    print("\t\t * merging temporary files\n");
    merge_RSCU(\%pars, \%data);

    print("\t\t* computing statistics and plotting results\n");
    plot_RSCU(\%pars, \%data);

  }

  elsif ($task =~/^DINUC_ODDS_RATIO$/i) {

    print("\t* computing dinucleotide odds ratio data\n");
    compute_dinuc(\%pars, \%data);

    print("\t\t* merging temporary files\n");
    merge_dinuc(\%pars, \%data);

    print("\t\t* computing statistics and plotting results\n");
    plot_dinuc(\%pars, \%data);

  }

  elsif ($task =~ /^ENC$/i) {
    
    print("\t* computing Effective Number of Codons (ENC) data\n");
    compute_ENC(\%pars, \%data);
  
    print("\t\t* merging temporary files\n");
    merge_ENC(\%pars, \%data);

    print("\t\t* computing statistics and plotting results\n");
    plot_ENC(\%pars, \%data);
  }
  
  elsif ($task =~ /^NUC_PER_POSITION$/i) {

    print("\t* computing Effective Number of Codons (ENC) data\n");
    compute_nuc_per_pos(\%pars, \%data);

    print("\t* computing Effective Number of Codons (ENC) data\n");
    merge_nuc_per_pos(\%pars, \%data);

    print("\t* computing Effective Number of Codons (ENC) data\n");
    plot_nuc_per_pos(\%pars, \%data);


  } else {
    die("Task $task not supported by SIGNATURES, please check if you entered the right information. Currently supported tasks are @supported_tasks\n");
  }
}

#will take genbank or fasta files as input and generate
#clean sequence data (
sub clean_sequence_data {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $out_dir_path = $pars{outdir_path}."tmp/parsed_files/";
  my $log_file_path = $pars{outdir_path}."logs/";
  my $sequence_length_cutoff = $pars{sequence_length_cutoff} || 100;
  my $filter_par = "all"; #must go to configuration file
  foreach my $key (keys %data) {
    my $outfile = $out_dir_path.$data{$key}{fancy_name}."_nt.fasta";
    if (-e $outfile) {
      print ("$key\_nt.fasta already exists\n");
      next;
    }
    print("Parsing file $key\n");
    my $file_path = $data{$key}{file_path};
#    my $outfile = $out_dir_path."$data{$key}{fanc}"
    if ($file_path =~ /.gb/) {
      my $command = "bin/extract_ORF.pl $file_path $out_dir_path $log_file_path $filter_par $sequence_length_cutoff";
      system($command);
#      print($command."\n");
#      my $a = <STDIN>;
    }
  }
  plot_clean_sequence_data_log();
}

sub plot_clean_sequence_data_log {

}

sub compute_nuc_per_pos {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $out_dir_path = $pars{outdir_path}."tmp/nuc_per_pos_files/";
  my $input_dir_path = $pars{outdir_path}."tmp/parsed_files/";
  my $log_file_path = $pars{outdir_path}."logs/";
  my $outfile = $out_dir_path."final_results_nuc_per_pos.txt";
  if (-e $outfile) {
    print ("final_results_nuc_per_pos.txt already exists\n");
  }
  else {
    foreach my $key (keys %data) {
      my $infile = "$input_dir_path$data{$key}{fancy_name}\_nt.fasta";
      if (-e $infile) {
        my $command = "perl bin/compute_nuc_per_pos.pl $infile $out_dir_path $data{$key}{group}";
        system $command;
        print("$command\n");
      }
    }
  }
}

sub merge_nuc_per_pos {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $results_dir = $pars{outdir_path}."tmp/nuc_per_pos_files/";
  opendir(DIR, $results_dir);
  my @files = readdir(DIR);
  closedir(DIR);
  my $outfile_path = $results_dir."final_results_nuc_per_pos.txt";
  if (-e $outfile_path) {#remove final file from previous analysis
    #    print("Here!!!\n\n");
    #my $a = <STDIN>;
    my $command = "rm $outfile_path";
    system $command;
    #    print ("There!\n");
    #my $a = <STDIN>;
  }
  open(OUT, ">$outfile_path");
  print OUT ("GENOME\tLABEL\tNUCLEOTIDE\tPOSITION\tVALUE\n");
  my $i = 0;
  foreach my $file (@files) {
    next if (($file eq ".") || ($file eq "..") || ($file !~ /_nuc_per_codon_pos.txt$/));
    open(IN, "$results_dir$file");
    my $header = <IN>;
    while (my $line = <IN>) {
      print OUT $line;
    }
    close IN;
  print("$file\t$i\n\n");
  $i++;
  }
  close OUT;
}

sub plot_nuc_per_pos {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $infile = $pars{outdir_path}."tmp/nuc_per_pos_files/final_results_nuc_per_pos.txt";
  my $outdir = $pars{outdir_path}."plots/";

  #computing main plot results
  my $command = "Rscript --vanilla --slave bin/nuc_per_pos_main_plot.R $infile $outdir"; #FALSE for plotting names (1st) and for Ncp plot (second)
  print($command."\n");
  system $command;

}

sub compute_ENC {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $out_dir_path = $pars{outdir_path}."tmp/ENC_files/";
  my $input_dir_path = $pars{outdir_path}."tmp/parsed_files/";
  my $log_file_path = $pars{outdir_path}."logs/";
  my $outfile = $out_dir_path."final_results_ENC.txt";
  if (-e $outfile) {
    print ("final_results_ENC.txt already exists\n");
  }
  else {
    foreach my $key (keys %data) {
      my $infile = "$input_dir_path$data{$key}{fancy_name}\_nt.fasta";
      #generating a concatenated CDS file 
      if (-e $infile) {
       #generating a concatenated CDS file
        my $command = "perl bin/concatenate_CDS.pl $infile $out_dir_path $key";
        system $command;
#        print $command."\t\t$key!!!\n";
#        my $a = <STDIN>;
        $infile = "$out_dir_path$data{$key}{fancy_name}\_nt_concatenated.fasta";
#        if ((defined $path2external{ENCprime})&&(defined $path2external{SeqCount})&&(-e $path2external{ENCprime})) {
          $command = "perl bin/compute_ENCprime.pl $path2external{ENCprime} $path2external{SeqCount} $infile $out_dir_path";
#          system $command;
          system $command."\n";
#          my $a = <STDIN>;
#        }
      }
      else {
    
      }
    }
  }
}

sub merge_ENC {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $results_dir = $pars{outdir_path}."tmp/ENC_files/";
  opendir(DIR, $results_dir);
  my @files = readdir(DIR);
  closedir(DIR);
  my $outfile_path = $results_dir."final_results_ENC.txt";
  open(OUT, ">$outfile_path");
  #  $header = "";
  my %tmp_data; #will store the parsed information about each genome
  print OUT ("GENOME\tLABEL\tNc\tNcp\tA_freq\tC_freq\tG_freq\tT_freq\n");
  foreach my $file (@files) {
    next if (($file eq ".") || ($file eq ".."));
    if ($file =~ /\.ENCprime\.txt$/) {
      open(IN, "$results_dir$file");
      my $header = <IN>;
      my $line = <IN>;
      chomp $line;
      my @aux = split(/\s+/, $line);
      my $genome = shift @aux;
      my $Nc = shift @aux;
      my $Ncp = shift @aux;
      $genome =~ s/:$//;
      $tmp_data{$genome}{Nc} = $Nc;
      $tmp_data{$genome}{Ncp} = $Ncp;
      close IN;
    }
    if ($file =~ /_nt_concatenated.fasta.acgtfreq$/) {
      open(IN, "$results_dir$file");
      my $header = <IN>;
      my $line = <IN>;
      chomp $line;
      my @aux = split(/\s+/, $line);
      my $genome = shift @aux;
      my $a = shift @aux;
      my $c = shift @aux;
      my $g = shift @aux;
      my $t = shift @aux;
      $genome =~ s/>$//;
#      print()
      $tmp_data{$genome}{a} = $a;
      $tmp_data{$genome}{c} = $c;
      $tmp_data{$genome}{g} = $g;
      $tmp_data{$genome}{t} = $t;
      close IN;
    }

  }
  foreach my $key (keys %tmp_data) {
    print OUT ("$key\t$data{$key}{group}\t$tmp_data{$key}{Nc}\t$tmp_data{$key}{Ncp}\t$tmp_data{$key}{a}\t$tmp_data{$key}{c}\t$tmp_data{$key}{g}\t$tmp_data{$key}{g}\n");
  }
  close OUT;
}

sub plot_ENC {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $infile = $pars{outdir_path}."tmp/ENC_files/final_results_ENC.txt";
  my $outdir = $pars{outdir_path}."plots/";

  #computing main plot results
  my $command = "Rscript --vanilla --slave bin/ENC_main_plot.R $infile $outdir FALSE FALSE"; #FALSE for plotting names (1st) and for Ncp plot (second)
  print($command."\n");
  system $command;

  $command = "Rscript --vanilla --slave bin/ENC_main_plot.R $infile $outdir FALSE TRUE"; #FALSE for plotting names (1st) and TRUE for Ncp plot (second)
  print($command."\n");
  system $command;

  $command = "Rscript --vanilla --slave bin/ENC_plot_by_group.R $infile $outdir FALSE FALSE"; #FALSE for plotting names (1st) and TRUE for Ncp plot (second)
  print($command."\n");
  system $command;

  $command = "Rscript --vanilla --slave bin/ENC_plot_by_group.R $infile $outdir FALSE TRUE"; #FALSE for plotting names (1st) and TRUE for Ncp plot (second)
  print($command."\n");
  system $command;

  #computing heatmap results
  #  $command = "Rscript --vanilla --slave bin/dinuc_heatmap_plot.R $infile $outdir";
  #  system $command;
}

sub compute_dinuc {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $out_dir_path = $pars{outdir_path}."tmp/dinuc_files/";
  my $input_dir_path = $pars{outdir_path}."tmp/parsed_files/";
  #  my @files = readdir($input_dir_path);
  my $log_file_path = $pars{outdir_path}."logs/";

  my $outfile = $pars{outdir_path}."tmp/dinuc_files/final_results_dinuc.txt";

  if (-e $outfile) {
    print ("final_results_dinuc.txt already exists\n");
  } else {
    foreach my $key (keys %data) {
      my $strandedness = $data{$key}{strandedness};
      my $strand_par;
      if ($strandedness =~ /^single$/i) {
        $strand_par = 0;
      } elsif ($strandedness =~ /^double$/i) {
        $strand_par = 1;
      } else {
        die ("Strandedness parameter should be either \"single\" or \"double\" for each genome, but you provided $strandedness for $key, please configure the metadata file according.\n");
      }
      my $infile = "$input_dir_path$data{$key}{fancy_name}\_nt.fasta";
#      print $infile."\n";
      if (-e $infile) {#if there is a parsed file
        my $command = "perl bin/dinucleotide_odds_ratio.pl $infile $out_dir_path $log_file_path $strand_par 0 $data{$key}{group}";
#        system($command);
        system $command;
#         my $command = "perl bin/dinucleotide_odds_ratio.pl $infile $out_dir_path $log_file_path $strand_par 0 $data{$key}{group} hclust";
#        my $a = <STDIN>;
      } else {
      }
    }
  }
}

sub merge_dinuc {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $results_dir = $pars{outdir_path}."tmp/dinuc_files/";
  opendir(DIR, $results_dir);
  my @files = readdir(DIR);
  closedir(DIR);
  my %tmp_file_data;
  my $outfile_path = $results_dir."final_results_dinuc.txt";
  open(OUT, ">$outfile_path");
  #  $header = "";
  print OUT ("GENOME\tLABEL\tDINUC\tOBSERVED\tEXPECTED\tODDS\n");
  foreach my $file (@files) {
    next if ($file !~ m/dinucleotide_odds_ratio.txt$/);
    open(IN, "$results_dir$file");
    my $header = <IN>;
    my @aux = split(/\t/, $header);
    if ($#aux == 6) {
      #      print("$header\n");
      while (my $line = <IN>) {
        chomp $line;
        my @aux = split(/\t/, $line);
        my $genome = $aux[0];
        my $label = $aux[1];
        my $dinuc = $aux[2];
        my $obs = $aux[4];
        my $exp = $aux[5];
        my $odds = $aux[6];
        print OUT ("$genome\t$label\t$dinuc\t$obs\t$exp\t$odds\n");
      }
    } elsif ($#aux == 8) {
      #      print("$header\n");
      while (my $line = <IN>) {
        chomp $line;
        my @aux = split(/\t/, $line);
        my $genome = $aux[0];
        my $label = $aux[1];
        my $dinuc = $aux[2];
        my $obs = $aux[6];
        my $exp = $aux[7];
        my $odds = $aux[8];
        print OUT ("$genome\t$label\t$dinuc\t$obs\t$exp\t$odds\n");
      }
    } else {
      die("merge_dinuc expects files with 7 or 9 columns, yours have $#aux, please delete your data and restart your analysis, your files may be broken\n");
    }
  }
  close OUT;
}

sub plot_dinuc {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $infile = $pars{outdir_path}."tmp/dinuc_files/final_results_dinuc.txt";
  my $outdir = $pars{outdir_path}."plots/";
  
  #computing main plot results
  my $command = "Rscript --vanilla --slave bin/dinuc_main_plot.R $infile $outdir";
  print($command."\n");
  system $command;
  
  #computing heatmap results
  $command = "Rscript --vanilla --slave bin/dinuc_heatmap_plot.R $infile $outdir";
  system $command;

  #computing pvclust results
  if ($pars{pvclust_boolean} eq "TRUE") {
    my $bootstrap = 1000;
    if (defined $pars{pvclust_bootstrap}) {
      $bootstrap = $pars{pvclust_bootstrap};
    }
    $command = "Rscript --vanilla --slave bin/dinuc_pvclust_plot.R $infile $outdir $bootstrap";
    print($command."\n");
    system $command;
  }

  #computing correlation matrix results, sort by alphabet
  $command = "Rscript --vanilla --slave bin/dinuc_corr_matrix_plot.R $infile $outdir alphabet";
  system $command;
  print($command."\n");
  #computing correlation matrix results, sort by hclust
  $command = "Rscript --vanilla --slave bin/dinuc_corr_matrix_plot.R $infile $outdir hclust";
  system $command;
  print($command."\n");

}

sub compute_RSCU {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $out_dir_path = $pars{outdir_path}."tmp/RSCU_files/";
  my $input_dir_path = $pars{outdir_path}."tmp/parsed_files/";
#    my @files = readdir($input_dir_path);
  my $log_file_path = $pars{outdir_path}."logs/";
  my $outfile = $pars{outdir_path}."tmp/RSCU_files/final_results_RSCU_minus_stop_monocodonic.txt";
  if (-e $outfile) {
    print ("final_results_RSCU_minus_stop_monocodonic.txt already exists\n");
  } else {
    foreach my $key (keys %data) {
      my $label = $data{$key}{group};
      my $infile = "$input_dir_path$data{$key}{fancy_name}\_nt.fasta";
      if (-e $infile) {#if there is a parsed file
        my $command = "perl bin/RSCU.pl $infile $out_dir_path $log_file_path $pars{codon_table} $label 0";
#        print $command."\n";
        system $command;
      } else {
      }
    }
  }
}
sub merge_RSCU {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $results_dir = $pars{outdir_path}."tmp/RSCU_files/";
  opendir(DIR, $results_dir);
  my @files = readdir(DIR);
  closedir(DIR);
  my $outfile_path = $results_dir."final_results_RSCU.txt";
  my $outfile2_path = $results_dir."final_results_RSCU_minus_stop_monocodonic.txt";
  open(OUT, ">$outfile_path");
  open(OUT2,">$outfile2_path");
  print OUT ("GENOME\tLABEL\tCODON\tAA\tRSCU\n");
  print OUT2 ("GENOME\tLABEL\tCODON\tAA\tRSCU\n");
  foreach my $file (@files) {
    next if ($file !~ m/_rel_syn_cod_usa.txt$/);
    open(IN, "$results_dir$file");
    my $header = <IN>;
    while (my $line = <IN>) {
      chomp $line;
      my @aux = split(/\t/, $line);
      my $genome = $aux[0];
      my $label = $aux[1];
      my $codon = $aux[2];
      my $degeneracy = $aux[6];
      my $aa = $aux[3];
      my $rscu = $aux[8];
#      print ("*$degeneracy*\t*$aa*\n");
#      my $a = <STDIN>;
      print OUT ("$genome\t$label\t$codon\t$aa\t$rscu\n");
      if (($degeneracy > 1)&&($aa ne "\*")) { #removing stop codons and monocodonic amino acids
        print OUT2 ("$genome\t$label\t$codon\t$aa\t$rscu\n");
      }
    }
  }
  close OUT;
  close OUT2;
}

sub plot_RSCU {
  my %pars = %{$_[0]};
  my %data = %{$_[1]};
  my $infile = $pars{outdir_path}."tmp/RSCU_files/final_results_RSCU_minus_stop_monocodonic.txt";
  my $outdir = $pars{outdir_path}."plots/";

  #computing heatmapt results
  my $command = "Rscript --vanilla --slave bin/RSCU_heatmap_plot.R $infile $outdir";
  print($command."\n");
  system $command;

  #computing pvclust
  if ($pars{pvclust_boolean} eq "TRUE") {
    my $bootstrap = 1000;
    if (defined $pars{pvclust_bootstrap}) {
      $bootstrap = $pars{pvclust_bootstrap};
    }
    my $command = "Rscript --vanilla --slave bin/RSCU_pvclust_plot.R $infile $outdir";
    print($command."\n");
    system $command;
  }

  #computing correlation matrix
  $command = "Rscript --vanilla --slave bin/RSCU_corr_matrix_plot.R $infile $outdir alphabet_aa";
  print($command."\n");
  system $command;

  #computing correlation matrix
  $command = "Rscript --vanilla --slave bin/RSCU_corr_matrix_plot.R $infile $outdir hclust";
  print($command."\n");
  system $command;

  #computing correlation matrix
  $command = "Rscript --vanilla --slave bin/RSCU_corr_matrix_plot.R $infile $outdir alphabet_nt";
  print($command."\n");
  system $command;
}

sub create_output_report {
  my  @tasts = @{$_[0]};
  my %pars = %{$_[1]};
  my %data = %{$_[2]};
}

sub validate_parameters {
  my %tmp_hash = %{$_[0]};
  if ($tmp_hash{outdir_path}) {
    my $outdir_path = $tmp_hash{outdir_path};
    if ((-d $outdir_path)) { #check if it's a writable directory
    } else {
      die("Please check if $outdir_path exists.\n");
    }
  } else {
    die("Please provide a path to the directory where SIGNATURES will write its output (\"outdir_path\" parameter in config_file).\n");
  }
  if ($tmp_hash{tasks}) {
    my $tasks = $tmp_hash{tasks};
  }
  else {
    die("Please provide a set of tasks for SIGNATURES to compute (\"tasks\" parameter in config_file).\n");
  }
}

#subroutines
sub usage {
  my $ver = $_[0];
  print "\n\n\t\tSIGNATURES version $ver.\n\n\n";
  print "Use this sofware like:\n";
  print "SIGNATURES <path to config_file.txt>\n";
  die();
}

sub parse_config_file {
  my $file_path = $_[0];
  my %tmp_data;
  if (-e $file_path) {
    open(IN, $file_path);
    while(my $line = <IN>) {
      next if (($line =~ /^#/)||($line =~ /^\n$/)||($line =~ /^\s+\n$/));#skipping comment and blank lines
      chomp $line;
      my @aux = split(/\s*=\s*/, $line);
      $tmp_data{$aux[0]} = $aux[1];
    }
    close IN;
    return \%tmp_data;
  }
  else {
    die ("File $file_path does not exist!\n");
  }
}

sub parse_metadata_file {
  my $tmp_file_path = $_[0];
  my %data;
  my $file_name = 0;
  if (-e $tmp_file_path) {
    open(IN, $tmp_file_path);
    my $header = "";
    while(my $line = <IN>) {
      chomp $line;
      next if (($line =~ /^#/)||($line =~ /^\n$/)||($line =~ /^\s+\n$/));#skipping comment and blank lines
      if ($line =~ /^GENOME\t/) {
        $header = $line;
      } else {
        my @aux = split(/\t/, $line);
        my @aux2 = split(/\//, $aux[0]);
        $file_name = pop @aux2;
        @aux2 = "";
        @aux2 = split(/\./, $file_name);
        my $fancy_name = $aux2[0]; #removing .gb or .fa from file name
        if (defined $data{$file_name}{fancy_name}) {
          die("There are another file named $fancy_name, please change it since SIGNATURE requires each genome/species to have a distinct name\n");
        }
        $data{$file_name}{fancy_name} = $fancy_name;
        #        $data{$file_name}
        $data{$file_name}{file_path} = $aux[0];
        $data{$file_name}{strandedness} = $aux[1];
        $data{$file_name}{group} = $aux[2];
        $data{$file_name}{status} = "OK"; #everything starts alright :-)
      }
    }
  } else {
    die("metadata file $tmp_file_path does not exist!\n");
  }
  return(\%data);
}

sub parse_path2external_file {
  my $file = $_[0];
  my %paths; #stores path to third-party software
  open(IN, $file) ||
    die("Could not open configuration file $file, please check configuration parameters");
  while (my $line = <IN>) {
    next if ($line =~ /^\n$/);
    next if ($line =~ /^\s*\n$/);
    next if ($line =~ /#/);
    chomp $line;
    my @aux = split(/\t/, $line);
    my $program = shift @aux;
    my $path = shift @aux;
    if ($program eq "ENCprime") {
      $paths{ENCprime} = $path;
    } elsif ($program eq "SeqCount") {
      $paths{SeqCount} = $path;
    } 
    else {
      print "$program does not match any third-software software currently supported by SIGNATURES, you may want to check what is happening here...\n";
    }
  }
  close IN;
  return(\%paths);
}
