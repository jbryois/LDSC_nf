#!/usr/bin/env perl

#Define hash of SNPs
my %snps_hash;

#File containing SNPs to be added to the LDscore annotations
my $snps=$ARGV[0];

#Name of the annotation
my $name=$ARGV[1];

#Reads file with SNPs to add to annotation and store them in a hash
open my $snps_file, '<', $snps or die "What have you done?! $!";
while (<$snps_file>) {
  chomp;                
  $snps_hash{$_} = 1;  
}
close $snps_file;

#Open Annotation file of LDSC
my $ld_scores_annot=$ARGV[2];
if ($ld_scores_annot =~ /.gz$/) {
open(IN, "gunzip -c $ld_scores_annot |") || die "can’t open pipe to $ld_scores_annot";
}
else {
open(IN, $ld_scores_annot) || die "can’t open $ld_scores_annot";
}

#Print header 
my $header =<IN>;
chomp $header;
print $header . "\t" . $name . "\n";

#if SNP in annotation is also in hash -> Add 1 in the new annotation, else add 0.
while (<IN>) {
  chomp;
  my ($chr,$bp, $snp,$cm,$base) = split '\t';
  my $line = $_;
  #print $line . "\n";
  if (exists($snps_hash{$snp})){
  	print $line . "\t" . "1" . "\n";
  }
  else {
  	print $line . "\t" . "0" . "\n";
  }	
}
close IN;