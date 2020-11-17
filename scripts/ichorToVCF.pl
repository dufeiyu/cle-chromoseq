#!/usr/bin/perl

use File::Basename;
use Getopt::Long;
use strict;

sub ratio2abundance;

my $refseq = "~/refdata/hg38/all_sequences.fa";

my $minCNAsize = 5000000;
my $minCNAabund = 10.0;
my $discardAbund = 5.0;
my $gender = '';

GetOptions("r=s" => \$refseq,
	   "g=s" => \$gender,
	   "s|minsize=i" => \$minCNAsize,
	   "f|minabund=f" => \$minCNAabund,
	   "m|discardabund=f" => \$discardAbund);

my $segsfile = $ARGV[0];

die "No gender specified!" if !$gender;

my $name = basename($segsfile,".segs.txt");

print <<EOF;
##fileformat=VCFv4.2
##source=ichorCNA
##FILTER=<ID=CNASize,Description="CNA size filter">
##FILTER=<ID=CNAAbundance,Description="CNA abundance filter"> 
##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=LOG2RATIO,Number=.,Type=Float,Description="Log2 ratio from CNV analysis">
##INFO=<ID=ABUNDANCE,Number=.,Type=Float,Description="Estimated abdunance from Log2 ratio from CNV analysis">
##INFO=<ID=CNBINS,Number=1,Type=Integer,Description="Number of CN bins in CNA call">
##INFO=<ID=POS,Number=1,Type=Integer,Description="Position of the variant described in this record">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">
##INFO=<ID=set,Number=1,Type=String,Description="Source VCF for the merged record in CombineVariants">
##INFO=<ID=CN,Number=.,Type=Integer,Description="Integer copy number estimate from CNV analysis">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FilterICHOR="minCNAsize:$minCNAsize; minCNAabund=$minCNAabund"
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$name
EOF

# now open CNV file and print as VCF records
open(CNV,"$segsfile") || die "cant open CNV file";
<CNV>; # remove header
while(<CNV>){
    chomp;

    my @F = split("\t",$_);

    # use the copy call to filter, unless Y, then use heuristic
    if ($gender eq 'male' and $F[1] eq 'chrY' and $F[5] < -0.38){
	$F[7] = 'HETD';
	$F[6]--;
    }
    
    next if $F[7] =~ /NEUT/;
    
    my $svtype = "DEL";
    $svtype = "DUP" if $F[5] > 0;
    my $pos = $F[2]-1;
    my $refnt = `samtools faidx $refseq $F[1]:$pos-$pos | tail -n 1`;
    chomp $refnt;
    my $svlen = ($F[3]-$F[2]+1);

    my $nl = 2;
    # if male and sex chromosome then adjust the normal and observed copy number
    if ($F[1] =~ /X|Y/ and $gender eq 'male'){
	$nl = 1;
	if ($F[5] < 0){
	    $F[6] = $nl-1;
	} else {
	    $F[6] = $nl+1;
	}
    }

    my $abund = sprintf("%.1f",ratio2abundance($F[6],$nl,$F[5]) * 100);

    # do not even report variants that are below $failAbund
    next if $abund < $discardAbund;
    
    my @filter = ();
    
    # Filter if too small or too low
    push @filter, 'CNASize' if $svlen < $minCNAsize;
    push @filter, 'CNAAbundance' if $abund < $minCNAabund;

    my $filter = join(',',@filter);
    $filter = 'PASS' if scalar @filter == 0;
    
    $svlen = -$svlen if $svtype eq 'DEL';
    print join("\t",$F[1],$F[2]-1,"ICHOR:$F[1]_$F[2]_$F[3]",$refnt,"<$svtype>",".",$filter,
		 join(";","SVTYPE=$svtype","LOG2RATIO=$F[5]","CN=$F[6]","ABUNDANCE=$abund","CNBINS=$F[4]","POS=".($F[2]-1),"END=$F[3]","SVLEN=$svlen","IMPRECISE"),
		 "GT","./."),"\n";
}

sub ratio2abundance {
    my ($cn,$nl,$l2r) = @_;
    my $t = ($nl * (2**$l2r - 1)) / ($cn - $nl);
    $t;
}
