#!/usr/bin/perl

# Generate exon-exon junction libraries with null model (by shifting) support


use warnings;
use strict;

use lib substr($0,0,rindex($0,'/')); 

require RETRODUP;

my $usage = $0." [-m null_mode -g genome_dir] annotation_file\n\n";
$usage   .= "Null modes examples: 1-1, 2-2, 100, -1000";

my $NO_ID  = "__NO_GENE_ID__";
my $EXTEND = 100;

my $gencode_file       = "";
my $mode               =  0; # predefined for true junction lib
my ($shift_s,$shift_e) = (0,0); # predefined for true junction lib
my $genome_dir         = "./";

my $n_arg = scalar(@ARGV); # number of parameters
for (my $i = 0;$i < $n_arg;$i++) { # get parameters of create_junction_library.pl
    my $arg = $ARGV[$i];
    if ($arg eq "-m") {
	if (++$i < $n_arg) { $mode       = $ARGV[$i]; }
    } elsif ($arg eq "-g") {
	if (++$i < $n_arg) { $genome_dir = $ARGV[$i]; } # dir updated
    } else {
	if (length($gencode_file) > 0) { # predefined par. above
	    print STDERR "Umbiguous input.\n";
	    print STDERR $usage,"\n";
	    exit;
	}
	$gencode_file = $arg;
    }
}

if ($mode =~ m/^(\d+)\-(\d+)$/) { # check mode par. format
    $shift_s =  $1;
    $shift_e = -$2;
} elsif ($mode =~ m/^(\-\d+)$/) {
    $shift_s = $shift_e = $1;
} elsif ($mode =~ m/^(\d+)$/) {
    $shift_s = $shift_e = $1;
} else {
    print STDERR "Unrecognized mode.\n";
    print STDERR $usage,"\n";
    exit;
}

if (length($gencode_file) == 0) { # check gencode_file par. format
    print STDERR "No annotation file is provided.\n";
    print STDERR $usage,"\n";
    exit;
}

my @genes              = ();
my %gene_exons         = (); # hash for gene_exons
my $FILE = RETRODUP::openFile($gencode_file);
if (!$FILE) { exit; }

while (my $line = <$FILE>) {
    if (substr($line,0,1) eq "#")     { next; }
    if (!($line =~ "protein_coding")) { next; }
    my @w = split(/\t/,$line);
    if ($w[2] eq "gene") {
	my $id = getId($w[8]);
	my @tmp = (getChromosome($w[0]));
	$gene_exons{$id} = \@tmp; # $id is the hash key. left-side now holds a reference to @tmp
	push(@genes,$id); # store geneID in @genes
    } elsif ($w[2] eq "exon") {
	my $id = getId($w[8]);
	my $arr = $gene_exons{$id}; # reference to @tmp, whose value is chromosome?
	if (!$arr) { next; }
	push(@$arr,$w[3] + $shift_s); # shift coordinates according to null model
	push(@$arr,$w[4] + $shift_e);
    }
}
close($FILE);

my @junctions      = ();
my %junction_names = ();

foreach my $gene (@genes) { # generate junctions and junction names?
    my $arr = $gene_exons{$gene}; # information: chr, shifted coordinates
    my $n = scalar(@$arr);
    if ($n <= 3) { next; } # if single exon, ignored
    if (int($n/2)*2 + 1 != $n) { next; }
    my $c = $$arr[0];

    my %tmp = ();
    for (my $i = 1;$i < $n;$i += 2) { $tmp{$$arr[$i]} = 1; }
    my @ss = sort { $a <=> $b } keys(%tmp);

    %tmp = ();
    for (my $i = 2;$i < $n;$i += 2) { $tmp{$$arr[$i]} = 1; }
    my @ee = sort { $a <=> $b } keys(%tmp);

    my ($ns,$ne) = (scalar(@ss),scalar(@ee));
    for (my $ie = 0;$ie < $ne;$ie++) {
	my $s = $ee[$ie];
	my $prefix = $c.":".($s - $EXTEND + 1)."-".$s."|0|"; 
	                         # EXTEND def: take length EXTEND in each of the two exons, then merge these two of length EXTEND.
	
	for (my $is = 0;$is < $ns;$is++) {
	    my $e = $ss[$is];
	    if ($e <= $s) { next; }
	    my $id = $prefix.$c.":".$e."-".($e + $EXTEND - 1);
	    if ($junction_names{$id}) {
		print STDERR "Duplicate junction detected ",$id," ...\n";
		$junction_names{$id} .= ":".$gene;
	    } else {
		push(@junctions,$id);
		$junction_names{$id} = $gene;
	    }
	}
    }
}

my ($seq,$chrom) = ("","");
foreach my $j (@junctions) { # output result into file
    my ($coor1,$ins,$coor2) = split(/\|/,$j);
    my ($c1,$s1,$e1) = split(/[\:\-]/,$coor1);
    my ($c2,$s2,$e2) = split(/[\:\-]/,$coor2);
    if ($c1 ne $chrom) {
	$chrom = $c1;
	$seq = parseChromosome($c1);
    }
    my $seq1 = substr($seq,$s1 - 1,$e1 - $s1 + 1);
    if ($c2 ne $chrom) {
	$chrom = $c2;
	$seq = parseChromosome($c2); # match to ref genome to extract sequence?
    }
    my $seq2 = substr($seq,$s2 - 1,$e2 - $s2 + 1);
    print ">",$junction_names{$j},"|splj|",$j,"|0\n";
    print $seq1,$seq2,"\n";
}

exit;

## subroutine definitions

sub parseChromosome # parse G1K reference
{
    my $chr = shift;
    my $file = $genome_dir."/chr".$chr;
    my $SEQFILE = RETRODUP::openFile($file.".fa.gz");
    if (!$SEQFILE) { $SEQFILE = RETRODUP::openFile($file.".fasta.gz"); }
    if (!$SEQFILE) { $SEQFILE = RETRODUP::openFile($file.".fa"); } # exist
    if (!$SEQFILE) { $SEQFILE = RETRODUP::openFile($file.".fasta"); }
    if (!$SEQFILE) { return ""; }
    my ($header,$seq) = ("","");
    while (my $line = <$SEQFILE>) {
        chomp($line);
        if (length($line) <= 0) { next; }
        if (substr($line,0,1) eq ">") {
            $header = substr($line,1);
            $seq = "";
        } else { $seq .= $line; }
    }
    close($SEQFILE);
    return $seq;
}

sub getId # return 'gene_id'_'gene_name'
{
    my $in = shift;
    my $ret = "";
    if ($in =~ m/gene_id\s\"(ENSG.+)\"\;\stranscript_id/) {
	$ret = $1;
    } else {
	return $NO_ID;
    }
    if ($in =~ m/gene_name\s\"(.+)\"\;\stranscript_type/) {
	$ret .= "_".$1;
    } else {
	return $NO_ID;
    }
    return $ret;
}

sub getChromosome # keep only the number, remove 'chr'
{
    my $in = shift;
    if (substr($in,0,3) eq "chr") { return substr($in,3); }
    return $in;
}
