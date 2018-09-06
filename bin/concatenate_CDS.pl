use strict;
use warnings;
use Bio::SeqIO;

my $infile = $ARGV[0];

my $outdir = $ARGV[1];

my $seq_name = $ARGV[2];

chomp $seq_name;

#print ("This is key!\t$seq_name\n\n");

#my $a = <STDIN>;

my $seqio_object = Bio::SeqIO->new(-file => $infile); 

my $superseq; #

while (my $seq_object = $seqio_object->next_seq) {
  my $seq = $seq_object->seq();
  $superseq = $superseq.$seq;
}

my $name = fancy_name($infile);

my $outfile = $outdir.$name."_concatenated.fasta";

#print("$outfile\n\n");

#my $a = <STDIN>;
#die();

open (OUT, ">$outfile");
print OUT ">$seq_name\n$superseq\n";
close OUT;

sub fancy_name {
  my $tmp = $_[0];
  my @aux = split(/\t/, $tmp);
  my @aux2 = split(/\//, $aux[0]);
  my $file_name = pop @aux2;
  $file_name =~ s/.fasta$//;
  return $file_name;
}
