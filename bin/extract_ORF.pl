#!/usr/bin/perl

################################################################################
##                                                                            ##
## Copyright 2018 UFMG                                                        ##
## Authors: Tarcisio Jose Domingos Coutinho/Francisco Pereira Lobo            ##
## this program is free software: you can redistribute it and/or modify       ##
## it under the terms of the GNU General Public License as published by the   ##
## Free Software Foundation, version 3 of the License.                        ##
##                                                                            ##
## extract_ORFs.pl is distributed in the hope that it will be useful,         ##
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

# Takes as input a genbank file and extracts valid fasta CDS files for
# downstream analyses.

use strict;
use warnings;
use Bio::SeqIO;

#$file_path $out_dir $log_file_path $filter_par $sequence_length_cutoff

my $file_path = shift;

my $output_dir = shift;

my $log_file = shift;

my $flags = shift;

my $len_cutoff = shift;

my %data;

if (! defined $len_cutoff) {
  $len_cutoff = 100;
}

if (! defined $flags) {
  $flags = "all";
}

#getting output file name and checking if it exists already
my $file = "";
$file = get_fancy_name($file_path);
my $outfile = "$output_dir/$file\_nt.fasta";

if (-e $outfile) {
  #  print("File $outfile already exists, please execute SIGNATURES with -f flag to overwrite temporary results\n");
  exit(1);
}

my $in = new Bio::SeqIO(-format => 'genbank',
                           -file => $file_path);
my $seq;

my $i = 1;
while ($seq = $in->next_seq()) {
  my $sequence = $seq->seq();
  my $start;
  my $end;
  my $organism;
  my $source;
  my $cds;
  my $flag = 0;
  for my $feat_object ($seq->get_all_SeqFeatures) {
    for my $tag ($feat_object->get_all_tags) {
      for my $value ($feat_object->get_tag_values($tag)) {
        if ($tag eq "pseudo") { #remove pseudogenes
          $flag = 1;
          next;
        }
      }
    if ( $feat_object->primary_tag eq 'CDS' ) {
      my $result = eval {
        my $cds_object = $feat_object->spliced_seq;
        $cds = $cds_object->seq;
      };
      unless ($result) {
#        print $@;
        next;
      }
    }
    if ($feat_object->primary_tag eq "source") {
      $source = $file;
	    $source =~s/\.gb.*$//g;
      $source =~ s/\s+/_/g;
      for my $tag ($feat_object->get_all_tags) {
        if ($tag eq "organism") {
          for my $value ($feat_object->get_tag_values($tag)) {
            $organism = $value;
          }
        }
      }
    $organism =~ s/\s+/_/g;
    $organism =~ s/\-/_/g;
    }
  }
    if ($feat_object->primary_tag eq "CDS") {
      next if check_cds($cds); #next if cds contains any problem
      my $id = "$source\_$i";
#      print OUT (">$source\_$i");
      for my $tag ($feat_object->get_all_tags) {
        if ($tag eq "protein_id") {
          for my $value ($feat_object->get_tag_values($tag)) {
            $data{$id}{protein_id} = $value;
#            print OUT ("||protein_id:$value");
#            $i++;
          }
        }
        if ($tag eq "organism") {
          for my $value ($feat_object->get_tag_values($tag)) {
            $data{$id}{organism} = $value;
#            print OUT ("||organism:$value");
          }
        }
        if ($tag eq "locus_tag") {
          for my $value ($feat_object->get_tag_values($tag)) {
            $data{$id}{locus_tag} = $value;
#            print OUT ("||locus_tag:$value");
          }
        }
        if ($tag eq "gene") {
          for my $value ($feat_object->get_tag_values($tag)) {
            $data{$id}{gene_id} = $value;
#            print OUT ("||gene:$value");
          }
        }
      }
      $data{$id}{cds} = $cds;
#      print OUT ("\n$cds\n");
    }
  }
  $i++;
}

my $entry_count = scalar keys %data; # number of valid sequences;

if ($entry_count >= 1) { #if at least one sequence is valid after parsing and filtering
  my $outfile = "$output_dir/$file\_nt.fasta";
  open (OUT, ">", "$outfile") ||
    die("Could not open file $outfile, please check!\n");
  foreach my $key (keys %data) {
    my $id = $key;
    print OUT (">$key");
    my $protein_id;
    my $locus_tag;
    my $gene_id;
    my $cds;
    my $flag = 0;
    if (defined $data{$key}{protein_id}) {
      print OUT ("|protein_id:$data{$key}{protein_id}");
    }
    if (defined $data{$key}{locus_tag}) {
      print OUT ("|locus_tag:$data{$key}{locus_tag}");
    }
    if (defined $data{$key}{gene_id}) {
      print OUT ("|gene_id:$data{$key}{gene_id}");
    }
    if (defined $data{$key}{cds}) {
      print OUT ("\n$data{$key}{cds}\n");
    } else {
      die("File $key has an empty sequence, but it (the sequence, of course) shoudn't have done so far! Something went really wrong here, so SIGNATURES can't keep going! :-(\n");
    }
  }
}

sub check_cds {
  my $tmp_seq = $_[0];
  if ($tmp_seq eq "") { #return error if empty
    return 1;
  }
  my $start_codon = substr($tmp_seq, 0, 3);
  my $stop_codon = substr($tmp_seq, -3);
  if ($start_codon !~ /ATG|GTG|CTG/i) {
#    print "Not valid start codon\t$tmp_seq\n";
    return (1);
  }
  if ($stop_codon !~ /TAA|TAG|TGA/i) {
#    print "Not valid stop codon\t$tmp_seq\n";
    return (1);
  }
  if (length($tmp_seq) % 3 != 0) {
#    print "length not multiple of three\t$tmp_seq\n";
    return (1);
  }
  if ($tmp_seq =~ /[^ACGT]/i) {
#    print "non-standard nucleotides\t$tmp_seq\n";
    return (1);
  }
  if (length($tmp_seq) < $len_cutoff) {
#    print "length smaller than cutoff\n";
    return (1);
  }
  return (0);
}

sub get_fancy_name {
  my $tmp = $_[0];
  my @aux = split(/\//, $tmp);
  $tmp = pop @aux;
  @aux = "";
  @aux = split(/\./, $tmp);
  $tmp = shift @aux;
  return $tmp;
}

