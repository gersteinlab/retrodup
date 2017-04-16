#!/bin/csh -f

# Extract unmapped reads and re-map them to the exon-exon junctions

set BWA_EXE      = "bwa"
set BWA_OPTIONS  = "aln -q 15 -n 0.04 -o 0 -t 8"
set SAMTOOLS_EXE = "samtools"
set GZIP_EXE     = "gzip"
set PERL_EXE     = "perl"
set XARGS_EXE    = "xargs"
set TMP_DIR      = "/tmp"


set change_path = "Please change the path to execulate in this script."
if (`which $BWA_EXE >& /dev/null; echo $?` == 1) then
  echo "Can't run bwa!"
  echo $change_path
  exit;
endif
if (`which $SAMTOOLS_EXE >& /dev/null; echo $?` == 1) then
  echo "Can't run samtools!"
  echo $change_path
  exit;
endif
if (`which $GZIP_EXE >& /dev/null; echo $?` == 1) then
  echo "Can't run gzip!"
  echo $change_path
  exit;
endif
if (`which $PERL_EXE >& /dev/null; echo $?` == 1) then
  echo "Can't run perl!"
  echo $change_path
  exit;
endif
if (`which $XARGS_EXE >& /dev/null; echo $?` == 1) then
  echo "Can't run xargs!"
  echo $change_path
  exit;
endif

# Usage
set usage = "Usage: \
   "$0" -o output_suffix \
        -i genome_index_prefix1 [-i genome_index_prefix2 -i ...] \
        file1.bam [file2.bam ...]"

# Input check
set n_args = $#argv
if ($n_args == 0) then
  echo $usage:q
  exit;
endif

# Parsing input
set genomes = ()
set files = ()
set i = 1
while ($i <= $n_args)
 set arg = $argv[$i]
 if ("$arg" == "-i" || "$arg" == "-o") then
@  i = $i + 1
   if ("$arg" == "-i" && $i <= $n_args) then
     set genomes = ($genomes $argv[$i])
   endif
   if ("$arg" == "-o" && $i <= $n_args) then
     set suff = $argv[$i]
   endif
 else
    set files = ($files $arg)
 endif
@ i = $i + 1
end

# Checking suffix
if ($?suff == 0) then
  echo "No output suffix provide!"
  echo $usage:q
  exit;
endif

# Checking genome indices
if ($#genomes <= 0) then
  echo "No genome indices provide!"
  echo $usage:q
  exit;
endif
foreach g ($genomes)
    if (`ls $g.* >& /dev/null; echo $?` == 1) then
	echo "No files resembling index "$g" found!"
	echo "Please check that provided prefix is correct."
	exit
    endif
end

# Checking input files
if ($#files <= 0) then
 echo "No input files provide!"
 echo $usage:q
 exit;
endif
foreach f ($files)
    if (! -e $f) then
      echo "File "$f" doesn't exist ..."
      exit
    else
      if (! -r $f) then
        echo "File "$f" isn't readable ..."
        exit
      endif
    endif
end

# Generating name for temporary file
set fq  = $TMP_DIR"/rdup_"`date '+%S'``date '+%N'`".fastq.gz"
if (-e $fq) then
  echo "File "$fq" exists ..."
  set fq = $TMP_DIR"/rdup_"`date '+%S'``date '+%N'`".fastq.gz"
endif
if (-e $fq) then
  echo "File "$fq" exists!"
  echo "Can't generate unique name for .fastq file."
  echo "Aborting ..."
  echo "Please restart the script!"
  exit
endif

# Selcting secondary reads
echo "Creating file "$fq" with unmapped reads from the following files: "
foreach f ($files)
  echo "  "$f
end
ls $files | $XARGS_EXE -i $SAMTOOLS_EXE view -f 0x4 {} | \
            $PERL_EXE -ane ' ($flag,$seq,$qual) = ($F[1],$F[9],$F[10]); \
                             if ($flag & 0x10) {                        \
                               $seq  = reverse($seq);                   \
                               $seq  =~ tr/[ACTGactg]/[TGACtgac]/;      \
                               $qual = reverse($qual);                  \
                             }                                          \
                             print "\@",$F[0],"\n";                     \
                             print $seq,"\n";                           \
                             print "+\n";                               \
                             print $qual,"\n";' | \
            $GZIP_EXE > $fq

# Remapping
set alis = ()
foreach g ($genomes)
    if (`ls $g.* >& /dev/null; echo $?` == 1) then
	echo "No files resembling index "$g" found!"
	echo "Please check that provided prefix is correct."
    else
	echo "Aligning to genome index "$g" ..."
        set prefix = `echo $g | rev | cut -f 1 -d "/" | rev`
        set out = $prefix"."$suff".sam.gz"
	$BWA_EXE $BWA_OPTIONS $g $fq | $BWA_EXE samse $g - $fq | \
	$SAMTOOLS_EXE view -q 20 -F 0x4 -S - | \
	    $GZIP_EXE > $out
        set alis = ($alis $out)
    endif
end
echo "Created the following files with alignments:"
foreach o ($out)
    echo "  "$o
end

# Deleting temporary file
echo "Deleting file "$fq" ..."
rm $fq

exit;
