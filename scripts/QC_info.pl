#! /usr/bin/perl

#Copyright (C) 2018 Feiyu Du <fdu@wustl.edu>
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

die "Provide chromoseq file, wildcards or directory" if @ARGV == 0;

my @files;
for my $arg (@ARGV) {
    if (-d $arg) {
        my @fs = glob("$arg/*.chromoseq.txt");
        push @files, @fs;
    }
    elsif (-f $arg) {
        die "$arg is not a valid file" unless -s $arg;
        push @files, $arg;
    }
    elsif ($arg =~ /\*/) {
        my @fs = glob($arg);
        unless (@fs > 0) {
            die "$arg doesn not have any file";
        }
        push @files, @fs;
    }
    else {
        die "Unknown $arg";
    }
}

#my $out_file = './QC_info.txt';
#my $out_fh = IO::File->new(">$out_file") or die "fail to write to $out_file";

#$out_fh->print(join "\t", "SAMPLE", "AVERAGE_COVERAGE", "TOTAL_READS", "MAPPED_READS", "MAPPED%", "DUPLICATES", "DUPLICATES%", "UNIQUE_READS", "UNIQUE%", "PROPERLY_PAIRED_READS", "PAIRED%", "MEAN_INSERT_SIZE", "CNAs", "RecurrentSVs", "GeneMutations", "OtherSVs");
#$out_fh->print("\n");

print join "\t", "SAMPLE", "AVERAGE_COVERAGE", "TOTAL_READS", "TOTAL_GIGABASES", "MAPPED_READS", "DUPLICATES", "UNIQUE_READS", "PROPERLY_PAIRED_READS", "MEAN_INSERT_SIZE", "CNAs", "RecurrentSVs", "GeneMutations", "OtherSVs";
print "\n";

for my $file (@files) {
    my $base = basename $file;
    my ($name) = $base =~ /^(\S+)\.chromoseq\.txt/;
    
    my ($av, $tr, $tg, $mr, $dup, $ur, $ppr, $mis);
    my ($cna_flag, $rt_flag, $gm_flag, $os_flag);
    my (@cna, @rt, @gm, @os);

    my $fh = IO::File->new($file) or die "failed to open $file for reading";
    while (my $line = $fh->getline) {
        if ($line =~ /\*\*\* COPY NUMBER/) {
            $cna_flag = 1;
        }
        elsif ($line =~ /\*\*\* RECURRENT TRANSLOCATIONS/) {
            $rt_flag  = 1;
            $cna_flag = 0;
        }
        elsif ($line =~ /\*\*\* GENE MUTATIONS/) {
            $gm_flag = 1;
            $rt_flag = 0;
        }
        elsif ($line =~ /\*\*\* FILTERED GENE MUTATIONS/) {
            $gm_flag = 0;
        }
        elsif ($line =~ /\*\*\* OTHER STRUCTURAL/) {
            $os_flag = 1;
        }
        elsif ($line =~ /\*\*\* CHROMOSEQ CASE INFORMATION/) {
            $os_flag = 0;
        }
        elsif ($line =~ /AVERAGE COVERAGE:\t(\S+)\t/) {
            $av = $1;
        }
        elsif ($line =~ /TOTAL READS:\t(\S+)\t/) {
            $tr = $1;
        }
        elsif ($line =~ /TOTAL GIGABASES:\t(\S+)\t/) {
            $tg = $1;
        }
        elsif ($line =~ /MAPPED READS:\t(\S+)\t/) {
            $mr = $1;
        }
        elsif ($line =~ /DUPLICATES:\t(\S+)\t/) {
            $dup = $1;
        }
        elsif ($line =~ /UNIQUE READS:\t(\S+)\t/) {
            $ur = $1;
        }
        elsif ($line =~ /PROPERLY PAIRED READS:\t(\S+)\t/) {
            $ppr = $1;
        }
        elsif ($line =~ /MEAN INSERT SIZE:\t(\S+)\t/) {
            $mis = $1;
        }

        if ($cna_flag) {
            my $cna = parse_sv($line);
            push @cna, $cna if $cna; 
        }
        if ($rt_flag) {
            my $rt = parse_sv($line);
            push @rt, $rt if $rt;
        }
        if ($gm_flag) {
            my $gm = parse_gm($line);
            push @gm, $gm if $gm;
        }
        if ($os_flag) {
            my $os = parse_sv($line);
            push @os, $os if $os;
        }
    }
    $fh->close;
    
    @cna = ('NA') unless @cna;
    @rt = ('NA') unless @rt;
    @gm = ('NA') unless @gm;
    @os = ('NA') unless @os;

    my $cna_str = join ',', @cna;
    my $rt_str = join ',', @rt;
    my $gm_str = join ',', @gm;
    my $os_str = join ',', @os;

    #$out_fh->print(join "\t", $name, $av, $tr, $mr, $mr_pct, $dup, $dup_pct, $ur, $ur_pct, $ppr, $ppr_pct, $mis, $cna_str, $rt_str, $gm_str, $os_str);
    #$out_fh->print("\n");
    print join "\t", $name, $av, $tr, $tg, $mr, $dup, $ur, $ppr, $mis, $cna_str, $rt_str, $gm_str, $os_str;
    print "\n";
}

#$out_fh->close;

sub adjust {
    my $num = shift;
    $num = 'NA' unless $num;
    my @nums = split /\s+/, $num;
    return @nums;
}

sub parse_sv {
    my $l = shift;
    return if $l =~ /^(TYPE|\*\*\*|\s+|\n)/;
    my @cols = split /\t/, $l;
    $cols[9] =~ s/seq\[GRCh38\]\s+//g;
    return $cols[9].'{'.$cols[7].'}['.$cols[13].']';
}

sub parse_gm {
    my $l = shift;
    return if $l =~ /^(TYPE|\*\*\*|\s+|\n)/;
    my @cols = split /\t/, $l;
    return $cols[5].':'.$cols[8].'['.$cols[12].']';
}
