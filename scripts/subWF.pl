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

use Cwd qw(abs_path);
use JSON qw(from_json to_json);
use IO::File;
use File::Spec;
use File::Basename;

die "Provide Chromoseq output directory" unless @ARGV == 1;
my $dir = $ARGV[0];
die "$dir is not a valid directory" unless -d $dir;
$dir = abs_path($dir);

my $git_dir = '/storage1/fs1/duncavagee/Active/SEQ/Chromoseq/process/git/cle-chromoseq';
my $conf = File::Spec->join($git_dir, 'application.conf');
my $wdl  = File::Spec->join($git_dir, 'ChromoseqAnalysis.wdl');
my $json_template = File::Spec->join($git_dir, 'inputs.analysis.json');

my $group  = '/cle/wdl/chromoseq';
my $queue  = 'dspencer';
my $docker = 'registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-20';

my $user_group = 'compute-gtac-mgi';

my $main_json = File::Spec->join($dir, 'test_inputs.json');  #Sometimes only part of cases need run subWF
die "$main_json is not valid" unless -s $main_json;

my $main_inputs = from_json(`cat $main_json`);
my $ct = scalar @{$main_inputs->{'ChromoSeq.Samples'}};
my %info;

for my $idx (0..$ct-1) {
    $info{$main_inputs->{'ChromoSeq.Samples'}->[$idx]} = {
        DOB       => $main_inputs->{'ChromoSeq.DOBs'}->[$idx],
        sex       => $main_inputs->{'ChromoSeq.Genders'}->[$idx],
        exception => $main_inputs->{'ChromoSeq.Exceptions'}->[$idx],
    };
}

for my $case_dir (glob("$dir/TW*")) {
    my ($case_name) = basename($case_dir) =~ /^(TW\S+)$/;
    unless ($info{$case_name}) {
        die "$case_name does not have info";
    }

    my $dragen_dir = File::Spec->join($case_dir, 'dragen');
    die "$dragen_dir not existing" unless -d $dragen_dir;
    
    my $inputs = from_json(`cat $json_template`);
    $inputs->{'ChromoseqAnalysis.Cram'}           = File::Spec->join($dragen_dir, $case_name.'_tumor.cram');
    $inputs->{'ChromoseqAnalysis.CramIndex'}      = File::Spec->join($dragen_dir, $case_name.'_tumor.cram.crai'); 
    $inputs->{'ChromoseqAnalysis.MappingSummary'} = File::Spec->join($dragen_dir, $case_name.'.mapping_metrics.csv');
    $inputs->{'ChromoseqAnalysis.CoverageSummary'}= File::Spec->join($dragen_dir, $case_name.'.wgs_coverage_metrics.csv');
    $inputs->{'ChromoseqAnalysis.TumorCounts'}    = File::Spec->join($dragen_dir, $case_name.'.target.counts.gz');

    $inputs->{'ChromoseqAnalysis.Name'}           = $case_name;
    $inputs->{'ChromoseqAnalysis.OutputDir'}      = $case_dir;

    $inputs->{'ChromoseqAnalysis.DOB'}            = $info{$case_name}->{DOB};
    $inputs->{'ChromoseqAnalysis.Gender'}         = $info{$case_name}->{sex};
    $inputs->{'ChromoseqAnalysis.Exception'}      = $info{$case_name}->{exception};
    $inputs->{'ChromoseqAnalysis.RunInfoString'}  = $main_inputs->{'ChromoSeq.RunInfoString'};

    my $input_json = File::Spec->join($case_dir, 'inputs.analysis.json');
    my $json_fh = IO::File->new(">$input_json") or die "fail to write to $input_json";

    $json_fh->print(to_json($inputs, {canonical => 1, pretty => 1}));
    $json_fh->close;

    my $out_log = File::Spec->join($case_dir, 'out.log');
    my $err_log = File::Spec->join($case_dir, 'err.log');

    my $cmd = "bsub -g $group -G $user_group -oo $out_log -eo $err_log -q $queue -R \"select[mem>8000] rusage[mem=8000]\" -M 8000000 -a \"docker($docker)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl -i $input_json $wdl";

    system $cmd;
    #print $cmd."\n";
    print $case_name." submitted\n";
}
