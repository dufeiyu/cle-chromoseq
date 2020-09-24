#!/usr/bin/perl

use strict;
use JSON;
use File::Spec;
use Getopt::Long;

my $base_dir = '/gscmnt/gc13016/cle/54f8f7b915cb472aa183c721307369ab_scratch_space/ChromoSeq/git/cle-chromoseq';
my $json = File::Spec->join($base_dir, 'inputs.local.json');
my $wdl  = File::Spec->join($base_dir, 'Chromoseq.wdl'); 
my $conf = File::Spec->join($base_dir, 'application.conf');

my $group  = '/cle/wdl/chromoseq';
my $queue  = 'compute-mgi-cle';
my $docker = 'registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:7';

my $dir = '';
my $name = '';
my $out = '';
my $gender = '';

my $testing = 0;

GetOptions(
    "gender=s" => \$gender,
    "i|inputs=s" => \$json,
    "w|wdl=s" => \$wdl,
    "d|dir=s" => \$dir,
    "n|name=s" => \$name,
    "o|out=s" => \$out,
    "g|group=s" => \$group,
    "t" => \$testing
);

die "$0 -d dir -n name -o out" if !-d $dir and !$name and !-d $out;
die "Input json, wdl, conf is(are) not valid" unless -s $json and -s $wdl and -s $conf;

my $inputs = from_json(`cat $json`);

chomp ($wdl = `readlink -f $wdl`);

chomp (my $mapsum = `readlink -f $dir/*.mapping_metrics.csv`);
die "Fail to find mapping_metrics.csv" if ! -e $mapsum;
$inputs->{'ChromoSeq.MappingSummary'} = $mapsum;

chomp (my $covsum = `readlink -f $dir/*.wgs_coverage_metrics.csv`);
$inputs->{'ChromoSeq.CoverageSummary'} = $covsum if -e $covsum;

chomp (my $cram = `readlink -f $dir/*.cram`);
die "Fail to find cram" if ! -e $cram;
$inputs->{'ChromoSeq.Cram'} = $cram;

chomp (my $crai = `readlink -f $cram.crai`);
die "Fail to find cram index crai" if ! -e $crai;
$inputs->{'ChromoSeq.CramIndex'} = $crai;

chomp (my $counts = `readlink -f $dir/*.target.counts`);
die "Fail to get tumor counts" if ! -e $counts;
$inputs->{'ChromoSeq.TumorCounts'} = $counts;

$inputs->{'ChromoSeq.Name'} = $name;

chomp ($out = `readlink -f $out`);
mkdir $out if ! -e $out;
$inputs->{'ChromoSeq.OutputDir'} = $out;

if ($gender =~ /^male$|^female$/){
    $inputs->{'ChromoSeq.Gender'} = $gender;

} 
elsif (`grep -c estimation $dir/*.ploidy_estimation_metrics.csv` =~ /1/) {
    `grep estimation $dir/*.ploidy_estimation_metrics.csv | cut -d ',' -f 4` =~ /(XX|XY)/;
    if ($1 eq 'XX'){
        $inputs->{'ChromoSeq.Gender'} = "female";
    } 
    elsif ($1 eq 'XY') {
        $inputs->{'ChromoSeq.Gender'} = "male";
    }
} 
elsif (`grep -c XX $dir/*.wgs_ploidy.csv` =~ /1/){
    $inputs->{'ChromoSeq.Gender'} = "female";

} 
elsif (`grep -c XY $dir/*.wgs_ploidy.csv` =~ /1/){
    $inputs->{'ChromoSeq.Gender'} = "male";
    
} 
elsif (`grep -c XAvgCov/YAvgCov $dir/*.wgs_coverage_metrics.csv` =~ /1/) {
    `grep XAvgCov/YAvgCov $dir/*.wgs_coverage_metrics.csv | cut -d ',' -f 4` =~ /([0-9.]+)/;
    if ($1 > 1.75){
        $inputs->{'ChromoSeq.Gender'} = "female";
    } 
    elsif ($1 < 1.75) {
        $inputs->{'ChromoSeq.Gender'} = "male";
    }
} 
elsif (`grep -c XAvgCov/YAvgCov $dir/*.mapping_metrics.csv` =~ /1/) {
    `grep XAvgCov/YAvgCov $dir/*.mapping_matrics.csv | cut -d ',' -f 4` =~ /([0-9.]+)/;
    if ($1 > 1.75){
        $inputs->{'ChromoSeq.Gender'} = "female";
    } 
    elsif ($1 < 1.75) {
        $inputs->{'ChromoSeq.Gender'} = "male";
    }
}

die "Cant find gender\n" if $inputs->{'ChromoSeq.Gender'} eq '';

open(JS,">$out/chromoseq.$name.json") || die;
print JS to_json($inputs);
close JS;

my $cmd = "bsub -g $group -oo $out/out.$name.log -eo $out/err.$name.log -q $queue -a \"docker($docker)\" /usr/bin/java -Dconfig.file=$conf -jar /opt/cromwell.jar run -t wdl -i $out/chromoseq.$name.json $wdl";

#print STDERR $cmd,"\n";

`$cmd` if !$testing;

print "$name WDL submitted\n";
