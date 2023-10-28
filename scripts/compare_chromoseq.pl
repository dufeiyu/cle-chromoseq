#! /usr/bin/perl

#Copyright (C) 2021 Feiyu Du <fdu@wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;
use IO::File;
use File::Basename;

umask 002;

die "Provide base and query chromoseq.txt files in order" unless @ARGV == 2;
my ($base_txt, $query_txt) = @ARGV;
die "$base_txt and or $query_txt is not valid" unless -s $base_txt and -s $query_txt;

my ($name) = basename($base_txt) =~ /^(\S+)\.chromoseq\.txt/;
die "$base_txt is not valid chromoseq.txt file" unless $name;

my @header  = ("Case");
my @compare = ($name);

my %base  = parse_txt($base_txt);
my %query = parse_txt($query_txt);

for my $category (sort keys %base) {
    push @header, $category.'_recall', $category.'_precision';
    my @base   = @{$base{$category}};
    my @query  = @{$query{$category}};
    if (@base == 0 and @query == 0) {
        push @compare, 'NA', 'NA';
    }
    else {
        my $TP = 0;
        for my $base_str (@base) {
            $TP++ if grep{$_ eq $base_str}@query;
        }
        my ($recall, $precision);
        if ($TP == 0) {
            $recall    = 0;
            $precision = 0;
        }
        else {
            $recall = sprintf("%.3f", $TP/@base);
            $precision = sprintf("%.3f", $TP/@query);
        }
        push @compare, $recall, $precision;
    }
}

print join "\t", @header;
print "\n";
print join "\t", @compare;
print "\n";

sub parse_txt {
    my $txt = shift;
    my $fh = IO::File->new($txt) or die "Fail to open $txt for reading";
    
    my ($cna, $knownsv, $var, $novelsv) = (0) x 4;
    my (@CNA, @KnownSV, @Var, @NovelSV);
    
    while (my $line = $fh->getline) {
        if ($line =~ /COPY NUMBER ALTERATIONS/) {
            $cna = 1;
        }
        elsif ($line =~ /RECURRENT TRANSLOCATIONS/) {
            $cna = 0;
            $knownsv = 1;
        }
        elsif ($line =~ /\*\s+GENE MUTATIONS/) {
            $knownsv = 0;
            $var = 1;
        }
        elsif ($line =~ /FILTERED GENE MUTATIONS/) {
            $var = 0;
        }
        elsif ($line =~ /OTHER STRUCTURAL VARIANTS/) {
            $novelsv = 1;
        }
        elsif ($line =~ /CHROMOSEQ CASE INFORMATION/) {
            $novelsv = 0;
            last;
        }
        
        next if $line =~ /^(\s+|\*|TYPE|ChromoSeq)/;
        my @cols = split /\t/, $line;
        my $str = join '_', map{$cols[$_]}(0..4);

        if ($cna) {
            push @CNA, $str;
        }
        if ($knownsv) {
            push @KnownSV, $str;
        }
        if ($var) {
            push @Var, $str;
        }
        if ($novelsv) {
            push @NovelSV, $str;
        }
    }
    my %hash = (
        'Var' => \@Var,
        'CNA' => \@CNA,
        'KnownSV' => \@KnownSV,
        'NovelSV' => \@NovelSV,
    );
    return %hash;
}

