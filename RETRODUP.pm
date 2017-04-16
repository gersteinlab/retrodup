package RETRODUP;

use warnings;
use strict;

my $NO_ID = "__NO_GENE_ID__";

sub getGeneId
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

sub openFile
{
    my $file_name = shift;
    unless (-e $file_name) {
	print STDERR "File ",$file_name," does not exist.\n";
	return 0;
    }
    my $ext = getExtension($file_name);
    my $ret = 0;
    if ($ext eq "gz" || $ext eq "gzip") {
	$ret = open(INFILE,"gzip -cd $file_name | ");
    } else {
	$ret = open(INFILE,$file_name);
    }
    my $line = <INFILE>;
    if ($ret) { return \*INFILE; }
    print STDERR "Can't open file '",$file_name,"'.\n";
    return 0;
}

sub getExtension
{
    my $file_name = shift;
    my @tmp = split(/\./,$file_name);
    my $n = scalar(@tmp);
    if ($n <= 1) { return ""; }
    return $tmp[$n - 1];
}

1;
