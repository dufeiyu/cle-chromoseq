#! /usr/bin/perl

#Copyright (C) 2021 Feiyu Du <fdu@wustl.edu>
#              and Washington University The Genome Institute

#This script is distributed in the hope that it will be useful, 
#but WITHOUT ANY WARRANTY or the implied warranty of 
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
#GNU General Public License for more details.


use strict;
use warnings;

umask 002;

use lib "/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/perl5/lib/perl5";
use Spreadsheet::Read;
use File::Copy::Recursive qw(dircopy);
use Data::Dumper;
use JSON qw(from_json to_json);
use IO::File;
use File::Spec;
use File::Basename;

##THIS LAUNCHER SCRIPT NEEDS TO BE RUN ON DRAGEN NODE compute1-dragen-2 TO BE ABLE TO COPY RUNDIR
##TO LOCAL STAGING DRIVE
die "Provide rundir, excel sample spreadsheet, and batch name in order" unless @ARGV == 3;

my ($rundir, $sample_sheet, $batch_name) = @ARGV;
die "$sample_sheet is not valid" unless -s $sample_sheet;

my $staging_rundir = '/staging/runs/Chromoseq/rundir';
my $dir = '/storage1/fs1/duncavagee/Active/SEQ/Chromoseq';

if ($rundir =~ /$staging_rundir/) {
    print "$rundir is ready\n";
}
else {
    die "Run on dragen node compute1-dragen-2 !" unless -d $staging_rundir;
    #chdir($staging_rundir) or die "cannot change: $!\n";
    die "$rundir is not valid" unless -d $rundir;
    $rundir =~ s/\/*$//;
    dircopy($rundir, $staging_rundir) or die $!;
    #my $rv = system "/bin/cp -a $rundir ./";
    #my $rv = system "/usr/bin/rsync -avW $rundir $staging_rundir";
    #unless ($rv == 0) {
    #    die "Failed to cp $rundir to $staging_rundir\n";
    #}
    print "Copying to $staging_rundir is DONE";
    $rundir = File::Spec->join($staging_rundir, basename($rundir));
    #chdir($dir) or die "cannot change: $!\n";
}

my $git_dir = File::Spec->join($dir, 'process', 'git', 'cle-chromoseq');

my $conf = File::Spec->join($git_dir, 'application.conf');
my $wdl  = File::Spec->join($git_dir, 'Chromoseq.wdl');
my $json_template = File::Spec->join($git_dir, 'inputs.json');

my $group  = '/cle/wdl/chromoseq';
my $queue  = 'pathology';
my $docker = 'registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-8';

my $user_group = 'compute-duncavagee';

my $out_dir = File::Spec->join($dir, 'batchdir', $batch_name);
unless (-d $out_dir) {
    unless (mkdir $out_dir) {
        die "Failed to make directory $out_dir";
    }
}

my $dragen_ss  = File::Spec->join($out_dir, 'sample_sheet.csv'); 
my $ss_fh = IO::File->new(">$dragen_ss") or die "fail to write to $dragen_ss";
$ss_fh->print("[Data]\n");
$ss_fh->print("Lane,Sample_ID,Sample_Name,Sample_Project,index,index2\n");

my $data = Spreadsheet::Read->new($sample_sheet);
my $sheet = $data->sheet(1);

my %info;

for my $row ($sheet->rows()) {
    next if $row->[0] =~ /Run|Lane/i;
    unless ($row->[0] =~ /\d+/) {
        die "Lane number is expected, Check sample sheet spreadsheet";
    }
    my ($lane, undef, $lib, $sex, $index1, $index2, $exception) = @$row;
    $exception = 'NONE' unless $exception;

    unless ($sex =~ /^(male|female)$/) {
        die "$lib sex has to be either male or female, not:$sex";
    }

    my ($clean_index1) = $index1 =~ /([ACGT]+)/;
    my ($clean_index2) = $index2 =~ /([ACGT]+)/;

    unless (length($clean_index1) =~ /^(8|10)$/ and length($clean_index2) =~ /^(8|10)$/) {
        die "Provided $lib index1:$index1 and or index2:$index2 do not have 8 or 10 bases";
    }

    my $fix_index2 = rev_comp($clean_index2);
    $ss_fh->print(join ',', $lane, $lib, $lib, '', $clean_index1, $fix_index2);
    $ss_fh->print("\n");

    $info{$lib} = {
        sex => $sex,
        exception => $exception,
    };
}
$ss_fh->close;

my $inputs = from_json(`cat $json_template`);
$inputs->{'ChromoSeq.RunDir'} = $rundir;
$inputs->{'ChromoSeq.SampleSheet'} = $dragen_ss;
$inputs->{'ChromoSeq.Batch'} = $batch_name;

my @samples;
my @genders;
my @exceptions;

for my $sample (sort keys %info) {
    push @samples, $sample;
    push @genders, $info{$sample}->{sex};
    push @exceptions, $info{$sample}->{exception};
}

$inputs->{'ChromoSeq.Samples'} = \@samples;
$inputs->{'ChromoSeq.Genders'} = \@genders;
$inputs->{'ChromoSeq.Exceptions'} = \@exceptions;

my $input_json = File::Spec->join($out_dir, 'inputs.json');
my $json_fh = IO::File->new(">$input_json") or die "fail to write to $input_json";

$json_fh->print(to_json($inputs));
$json_fh->close;

my $out_log = File::Spec->join($out_dir, 'out.log');
my $err_log = File::Spec->join($out_dir, 'err.log');

my $cmd = "bsub -g $group -G $user_group -oo $out_log -eo $err_log -q $queue -a \"docker($docker)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl -i $input_json $wdl";

system $cmd;


sub rev_comp {
    my $index = shift;
    my $revcomp = reverse $index;
    $revcomp =~ tr/ATGCatgc/TACGtacg/;

    return $revcomp;
}
