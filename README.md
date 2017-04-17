# Retrodup

An intergated retroduplication caller based on 1) exon-exon junction and 2) discordant reads
Please refer to /*Landscape and Variation of Novel Retroduplications in 26 Human Populations*/. 


1. Create junction libraries
`create_junction_library.pl -m null_mode -g genome_dir annotation_file (gencode)`
2. Extract unmapped reads and re-map to junctions
`run_bwa.sh -o output_suffix junction_libs bam_files >& out.txt`
3. Call retroduplications from exon-exon junction
`call_retroduplications.pl annotation_file sam_file > result.txt`
4. Cluster the discordant reads to discover insertion sites	
`samtools view bamfile gene_coor | java Cluster gene_coor2 > insertionpoints.txt`
