#!/usr/bin/perl

# Call retroduplicatins from the exon-exon junction

use strict;
use warnings;

use lib substr($0,0,rindex($0,'/'));
require RETRODUP;

my $BWA     = "-bwa";
my $BOWTIE  = "-bowtie";
my $aligner = $BWA;
my $MIN_INTRON = 70;

my ($gencode,$sig_file) = ("","");
my @back_files = ();
foreach my $arg (@ARGV) {
    if ($arg eq $BOWTIE) { $aligner = $BOWTIE; }
    elsif ($arg eq $BWA) { $aligner = $BWA; }
    else {
	if (length($gencode) <= 0)     { $gencode  = $arg; }
	elsif (length($sig_file) <= 0) { $sig_file = $arg; }
	else                           { push(@back_files,$arg); }
    }
}

my %gene_hash     = ();
my $IN_GENCODE = RETRODUP::openFile($gencode);
if (!$IN_GENCODE) { exit; }
while (my $line = <$IN_GENCODE>) {
    if (substr($line,0,1) eq "#") { next; }
    my @w = split(/\t/,$line);
    if ($w[2] ne "gene") { next; }
    my $id = RETRODUP::getGeneId($w[8]);
    $gene_hash{$id} = $w[0].":".$w[3]."-".$w[4];
}
close($IN_GENCODE);

my $signal_hash = parseEvidence($sig_file);
my @back_hashs = ();
foreach my $f (@back_files) { push(@back_hashs,parseEvidence($f)); }

my ($best_d,$best_e,$best_n) = (1e+10,1e+10,-1);

for (my $d = 1;$d <= 15;$d++) { print "\t",$d; }
print "\n";
for (my $err = 0.005;$err < 0.051;$err += 0.005) {
    print $err;
    for (my $d = 1;$d <= 15;$d++) {
	my $n_sig = callRetrodups($signal_hash,$err,$d,0);
	my $n_back = 0;
	foreach my $h (@back_hashs) {
	    my $n_tmp = callRetrodups($h,$err,$d,0);
	    if ($n_tmp > $n_back) { $n_back = $n_tmp; }
	}
	if ($n_sig > 0) {
		printf("\t%d (%2.0f)",$n_sig,$n_back/$n_sig*100);
	}
	if ($n_back > 0) { next; }
	if ($n_sig > $best_n ||
	    ($n_sig == $best_n && $err <= $best_e && $d >= $best_d)) {
	    $best_n = $n_sig;
	    $best_e = $err;
	    $best_d = $d;
	}
    }
    print "\n";
}

print "Calling with ",$best_e," and ",$best_d,"\n";
callRetrodups($signal_hash,$best_e,$best_d,2);

exit;


sub callRetrodups
{
    my ($evidence_hash,$err,$d,$out) = @_;
    my $delta = 5;
    my %gene_junctions = ();
    foreach my $g (keys(%$evidence_hash)) {
	my $arr = $$evidence_hash{$g};
	my $n = scalar(@$arr);
	for (my $i = 0;$i < $n;$i += $delta) {
	    if ($$arr[$i] >= $d && $$arr[$i + 1] >= $d &&
		$$arr[$i + 2] <= $err) {
		my $hash = $gene_junctions{$g};
		if (!$hash) {
		    my %tmp = ();
		    $gene_junctions{$g} = $hash = \%tmp;
		}
		my $junc = $$arr[$i + 3];
		if (!$$hash{$junc}) { $$hash{$junc} = 1; }
		else                { $$hash{$junc}++; }
	    }
	}
    }
    my $n_calls = 0;
    foreach my $g (keys(%gene_junctions)) {
	my $j_hash = $gene_junctions{$g};
	my @jj     = keys(%$j_hash);
	my $n = scalar(@jj);
	my $n_sequential = 1;
	for (my $i1 = 0;$i1 < $n;$i1++) {
	    my ($c1,$s1,$e1) = split(/[\:\-]/,$jj[$i1]);
	    for (my $i2 = $i1 + 1;$i2 < $n;$i2++) {
		my ($c2,$s2,$e2) = split(/[\:\-]/,$jj[$i2]);
		if (overlap($s1,$e1,$s2,$e2) > 0) { next; }
		if ($$j_hash{$jj[$i1]} > 1 || $$j_hash{$jj[$i2]} > 1) {
		    $n_sequential = 2;
		}
		for (my $i3 = $i2 + 1;$i3 < $n;$i3++) {
		    my ($c3,$s3,$e3) = split(/[\:\-]/,$jj[$i3]);
		    if (overlap($s1,$e1,$s3,$e3) > 0 ||
			overlap($s2,$e2,$s3,$e3) > 0) { next; }
		    if ($$j_hash{$jj[$i1]} > 1 ||
			$$j_hash{$jj[$i2]} > 1 ||
			$$j_hash{$jj[$i3]} > 1) {
			$n_sequential = 3;
			last;
		    }
		}
		if ($n_sequential == 3) { last; }
	    }
	    if ($n_sequential == 3) { last; }
	}
	if ($n_sequential > 1) {
	    $n_calls++;
	    my $n_j = scalar(keys(%$j_hash));
	    my $n_r = 0;
	    foreach my $key (keys(%$j_hash)) { $n_r += $$j_hash{$key}; }
	    if ($out > 0) { print $g,"\t",$n_j,"\t",$n_r,"\n"; }
	    if ($out > 1) {            # if $out==0, not print out
		foreach my $key (keys(%$j_hash)) {
		    print "\t",$key,"\t",$$j_hash{$key},"\n";
		}
		my $arr = $$evidence_hash{$g};
		my $n = scalar(@$arr);
		for (my $i = 0;$i < $n;$i += $delta) {
		    for (my $add = 0;$add < $delta;$add++) {
			print "\t",$$arr[$i + $add];
		    }
		    print "\n";
		}
	    }
	}
    }
    return $n_calls;
}

sub parseEvidence
{
    my $file = shift;
    my %evidence_hash = ();

    my $IN = \*STDIN;
    if ($file ne "-") { $IN = RETRODUP::openFile($file); }
    if (!$IN) { return 0; }
    while (my $line = <$IN>) {
	chomp($line);
	my @fields = split(/\t/,$line);
	my ($len,$n_mis) = (0,0);
	if ($aligner eq $BOWTIE) {
	    $len = length($fields[4]);
	    if (defined($fields[7])) {
		my @tmp = split(/\,/,$fields[7]);
		$n_mis = scalar(@tmp);
	    }
	} elsif ($aligner eq $BWA) {
	    if ($fields[5] =~ /(\d+)M/) { $len = $1; }
	    foreach my $f (@fields) {
		if (substr($f,0,5) ne "MD:Z:") { next; }
		foreach my $c (split(//,substr($f,5))) {
		    if ($c =~ m/[A-z]/) { $n_mis++; }
		}
		last;
	    }
	}
	my @tmp    = split(/\|/,$fields[2]);
	my @genes  = split(/\:/,$tmp[0]);
	my ($c1,$s1,$e1) = split(/[\:\-]/,$tmp[2]);
	my $ins = $tmp[3];
	my ($c2,$s2,$e2) = split(/[\:\-]/,$tmp[4]);
	if ($c1 ne $c2) { next; }
	if ($s2 - $e1 - 1 < $MIN_INTRON) { next; }
	my $js = $e1 - $s1 + 1;
	my $je = $js + 1 + $ins;
	my ($s,$e) = ($fields[3],$fields[3]);
	if ($aligner eq $BOWTIE) { $s++; }
	$e = $s + $len - 1;
	my $d1 = $js - $s + 1;
	my $d2 = $e - $je + 1;
	if ($d1 <= 0 || $d2 <= 0) { next; }
	foreach my $g (@genes) {
	    my $arr = $evidence_hash{$g};
	    if (!$arr) {
		my @new_arr = ();
		$evidence_hash{$g} = $arr = \@new_arr;
	    }
	    push(@$arr,$d1);
	    push(@$arr,$d2);
	    push(@$arr,$n_mis/$len);
	    push(@$arr,$c1.":".$e1."-".$s2);
	    push(@$arr,$line);

	}
    }
    close($IN);

    return \%evidence_hash;
}

sub overlap
{
    my ($s1,$e1,$s2,$e2) = @_;
    my $s = $s1; if ($s2 > $s) { $s = $s2; }
    my $e = $e1; if ($e2 < $e) { $e = $e2; }
    return $e - $s;
}
