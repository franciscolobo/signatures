use strict;
use warnings;

my $ENCprime = $ARGV[0];

my $ENCprime_path = get_path($ENCprime);

my $SeqCount = $ARGV[1];

my $SeqCount_path = get_path($SeqCount);

my $infile = $ARGV[2];

my $outdir = $ARGV[3];

chomp $outdir;

my $outfile_root_filename = fancy_name($infile);

my $outfile_ENCprime = $outdir.$outfile_root_filename.".ENCprime.txt";

#print("$outfile_codcnt\t$outfile_ENCprime\t$outfile_conf\n");

#my $a = <STDIN>;

if (-e $outfile_ENCprime) {
 
} else {
  my $command = "cd $SeqCount_path; ./SeqCount -c $infile 1";
  system $command;
#  system $command;
  $command = "cd $SeqCount_path; ./Seqcount -n $infile 1";
  system$command;
  $command = "";
  $command = "cd $ENCprime_path; ./ENCprime $infile.codcnt $infile.acgtfreq 1 tmp 0 -q; mv tmp $outfile_ENCprime\n";
 system $command;
}

sub fancy_name {
  my $name = $_[0];  
  my @aux = split(/\//, $name);
  my $file_name = pop @aux;
  @aux = "";
  @aux = split(/\./, $file_name);
  my $fancy_name = $aux[0]; #removing .gb or .fa from file name
  $fancy_name =~ s/_nt_nt_concatenated$//;
  return($fancy_name);
}

sub get_path {
  my $tmp = $_[0];
  my @aux = split(/\//, $tmp);
  my $last = pop @aux;
  my $path = join("/", @aux);
  return($path);
}
