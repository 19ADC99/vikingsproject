#! /usr/bin/env perl

=head1 Description

  Reads as par1 a proteome (multiple fasta), as par 2 a list (txt file) 
  with entries that will be kept from the dataset.
  
  Result to STDOUT


=cut

#------------------------------------------------------------

sub usage( $ )
  {
    print "$_[0]\n";
    system("pod2text $0");
    exit(1);
  }

#------------------------------------------------------------

(-e "$ARGV[0]" && -e "$ARGV[1]") || &usage("input file\n");
print STDERR "* Reading $ARGV[0]\n";
%prot=&fasta2hash("$ARGV[0]");
print STDERR "* Selecting seqs...\n";
open (LIST, "$ARGV[1]");
@LIST=<LIST>;
close (LIST);
chomp(@LIST);

foreach $keep (@LIST)
 {
# $keep =~ s/[\+-]$//g;
  $keep =~ s/^>//g;
# print STDERR "List $keep\n";
 if (exists $prot{$keep})
  {
  print STDOUT ">$keep\n$prot{$keep}\n";
  $n++;
  }
  else
  {
  #print STDERR "$keep not found\n";
  }
 }
 
print STDERR "* $n saved (".(scalar(@LIST)).") in list.\n";
 
# SUB

#------------------------------------------------------------
# Reads in an entire (multiple) fasta file and returns a hash in which
# the keys are the identifiers of the sequences (without the '>')
# and the values are the sequences themselves. Created by cesim 14/01/2002
#

sub fasta2hash ( $ )
 {
  my ($file,$key,$value);
  my (%fasta_hash);
  $file=$_[0];
  open (IN,$file);
  while (<IN>)
   {
    chomp;
    if (/^>(\S+)/)
     {
      $key=$1;
 #     print STDERR "FH key $key\t"; #debug
     } #if (/^>(\w)$/)
    else
     {
      (defined $key) || die "File $file is not a fasta file!";
      s/\s+//g;
      $fasta_hash{$key}.=uc($_);
     } #else
   } #while (<IN>)
  close IN;
  return (%fasta_hash);
 } #fasta2hash ( $ )
