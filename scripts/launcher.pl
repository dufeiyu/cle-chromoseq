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

die "Provide rundir, excel sample spreadsheet, and batch name in order" unless @ARGV == 3;

my ($rundir, $sample_sheet, $batch_name) = @ARGV;
die "$rundir is not valid" unless -d $rundir;
die "$sample_sheet is not valid" unless -s $sample_sheet;

my $dir = '/storage1/fs1/duncavagee/Active/SEQ/Chromoseq';
my $git_dir = File::Spec->join($dir, 'process', 'git', 'cle-chromoseq');

my $conf = File::Spec->join($git_dir, 'application.conf');
my $wdl  = File::Spec->join($git_dir, 'Chromoseq.wdl');
my $zip  = File::Spec->join($git_dir, 'imports.zip');
my $json_template = File::Spec->join($git_dir, 'inputs.json');

my $group  = '/cle/wdl/chromoseq';
my $queue  = 'dspencer';
my $docker = 'mgibio/genome_perl_environment:compute1-38';

my $user_group = 'compute-duncavagee';

my $out_dir = File::Spec->join($dir, 'batchdir', $batch_name);
unless (-d $out_dir) {
    unless (mkdir $out_dir) {
        die "Failed to make directory $out_dir";
    }
}

#parsing CoPath dump for MRN, ACCESSION, DOB and Gender
my $copath_dump = '/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/daily_accession/WML_Daily_Accession_Log_CLE.csv';
my $dump_fh = IO::File->new($copath_dump) or die "fail to open $copath_dump for reading";
my %hash;

while (my $l = $dump_fh->getline) {
    next if $l =~ /^(Accession|,,,,)/;
    next unless $l =~ /ChromoSeq\s/;
    
    my @columns = split /,/, $l;
    my $id  = $columns[0];
    my $sex = $columns[4];

    if ($sex eq 'M') {
        $sex = 'male';
    }
    elsif ($sex eq 'F') {
        $sex = 'female';
    }
    else {
        #Gender U
        warn "WARN: unknown gender $sex found for $id";
    }

    $hash{$id} = {
        MRN => $columns[2],
        DOB => $columns[5],
        sex => $sex,
    };
}
$dump_fh->close;

#parse sample spreadsheet
my $data = Spreadsheet::Read->new($sample_sheet);
my $sheet = $data->sheet(1);

my %info;
my $str;

for my $row ($sheet->rows()) {
    next if $row->[0] =~ /Run|Lane/i;
    unless ($row->[0] =~ /\d+/) {
        die "Lane number is expected, Check sample sheet spreadsheet";
    }
    my ($lane_str, undef, $name, $index, $exception) = @$row;

    $name =~ s/\s+//g;
    my ($accession) = $name =~ /^(W\d+\-\d+)\-/;
    unless ($accession) {
        die "$name must include accession id";
    }
    unless (exists $hash{$accession}) {
        die "Fail to get info from CoPath daily Active log for $accession of $name";
    }

    my @lib_name = ('TWDY', $hash{$accession}->{MRN}, $name);
    unless ($name =~ /lib\d+$/) {
        push @lib_name, 'lib1';
    }
    my $lib = join '-', @lib_name;

    $lane_str =~ s/\s+//g;
    my @lanes = split /,/, $lane_str;

    my ($index1, $index2) = $index =~ /([ATGC]{10})\-([ATGC]{10})/;
    
    for my $lane (@lanes) {
        $str .= join ',', $lane, $lib, $lib, '', $index1, $index2;
        $str .= "\n";
    }

    $exception = 'NONE' unless $exception;
    $info{$lib} = {
        DOB => $hash{$accession}->{DOB},
        sex => $hash{$accession}->{sex},
        exception => $exception,
    };
}

## Get RunInfoString
my $run_xml = File::Spec->join($rundir, 'RunParameters.xml');
unless (-s $run_xml) {
    die "RunParameters.xml $run_xml is not valid";
}
my $xml_fh = IO::File->new($run_xml) or die "Fail to open $run_xml";
my ($runid, $R1cycle, $R2cycle, $index1cycle, $index2cycle, $fcmode, $wftype, $instr, $side);
my $fctype = 0;

while (my $line = $xml_fh->getline) {
    if ($line =~ /<RunId>(\S+)<\/RunId>/) {
        $runid = $1;
    }
    elsif ($line =~ / ReadName="Read1" Cycles="(\d+)" /) {
        $R1cycle = $1;
    }
    elsif ($line =~ / ReadName="Read2" Cycles="(\d+)" /) {
        $R2cycle = $1;
    }
    elsif ($line =~ / ReadName="Index1" Cycles="(\d+)" /) {
        $index1cycle = $1;
    }
    elsif ($line =~ / ReadName="Index2" Cycles="(\d+)" /) {
        $index2cycle = $1;
    }
    elsif ($line =~ /<InstrumentType>(\S+)<\/InstrumentType>/) {
        $wftype = $1;
    }
    elsif ($line =~ /<InstrumentSerialNumber>(\S+)<\/InstrumentSerialNumber>/) {
        $instr = $1;
    }
    elsif ($line =~ /<Side>(\S+)<\/Side>/) {
        $side = $1;
    }
    elsif ($line =~ /<Type>FlowCell<\/Type>/) {
        $fctype = 1;
    }
    if ($fctype) {
        if ($line =~ /<Mode>(\S+)<\/Mode>/) {
            $fcmode = $1;
            $fctype = 0;
        }
    }
}
$xml_fh->close;
my $run_info_str = join ',', $runid, $instr, $side, $fcmode, $wftype, $R1cycle, $index1cycle, $index2cycle, $R2cycle; 

## DRAGEN sample sheet
my $dragen_ss  = File::Spec->join($out_dir, 'sample_sheet.csv'); 
my $ss_fh = IO::File->new(">$dragen_ss") or die "fail to write to $dragen_ss";

if ($index1cycle == 19) {
    print "Index1 has 19 cycles!\n";
    $ss_fh->print("[Settings]\n");
    $ss_fh->print("OverrideCycles,Y151;I10N9;I10;Y151\n");
}

$ss_fh->print("[Data]\n");
$ss_fh->print("Lane,Sample_ID,Sample_Name,Sample_Project,index,index2\n");
$ss_fh->print($str);
$ss_fh->close;

## Input JSON
my $inputs = from_json(`cat $json_template`);

$inputs->{'ChromoSeq.Batch'}         = $batch_name;
$inputs->{'ChromoSeq.RunDir'}        = $rundir;
$inputs->{'ChromoSeq.SampleSheet'}   = $dragen_ss;
$inputs->{'ChromoSeq.RunInfoString'} = $run_info_str;

my @samples = (sort keys %info);

$inputs->{'ChromoSeq.Samples'}    = \@samples;
$inputs->{'ChromoSeq.DOBs'}       = [map{$info{$_}->{DOB}}@samples];
$inputs->{'ChromoSeq.Genders'}    = [map{$info{$_}->{sex}}@samples];
$inputs->{'ChromoSeq.Exceptions'} = [map{$info{$_}->{exception}}@samples];

my $input_json = File::Spec->join($out_dir, 'inputs.json');
my $json_fh = IO::File->new(">$input_json") or die "fail to write to $input_json";

$json_fh->print(to_json($inputs, {canonical => 1, pretty => 1}));
$json_fh->close;

my $out_log = File::Spec->join($out_dir, 'out.log');
my $err_log = File::Spec->join($out_dir, 'err.log');

my $cmd = "bsub -g $group -G $user_group -oo $out_log -eo $err_log -q $queue -R \"select[mem>16000] rusage[mem=16000]\" -M 16000000 -a \"docker($docker)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl --imports $zip -i $input_json $wdl";

system $cmd;
#print $cmd."\n";

