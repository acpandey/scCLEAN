Generate ccs reads from subreads file + downsample

`ccs subreads.bam ccs.bam`
`samtools view -s $DOWNSAMPLE_RATIO > ccs.downsampled.bam`

Remove primers

`lima --isoseq --dump-clips ccs.downsampled.bam primers.fasta output.bam`

10X 3' single cell GEX library uses 12bp UMI and 16bp barcode

`isoseq tag output.5p--3p.bam flt.bam --design T-12U-16B`
`isoseq3 refine flt.bam primers.fasta fltnc.bam --require-polya`

Isoseq correct + group dedup with the default parameters

`isoseq3 correct --barcodes 3M-february-2018.txt fltnc.bam fltnc.corrected.bam`
`samtools sort -CB fltnc.corrected.bam > fltnc.sorted.corrected.bam`
`isoseq3 groupdedup fltnc.sorted.correccted.bam dedup.bam`

Align with pbmm2

`pbmm2 align --preset ISOSEQ --sort dedup.bam GCF_000001405.39_GRCh38.p13_genomic.fna mapped.bam`

Collapse into unique isoforms

`isoseq3 collapse maped.bam collapsed.gff`

Sort output GFF, classify and refine results with Pigeon filter

`pigeon sort collapsed.gff -o sorted.gff`

Index both gtf and fasta file with `pigeon index` and `samtools faidx`

`pigeon classify -fl sorted.gff GCF_000001405.39_GRCh38.p13_genomic.gtf GCF_000001405.39_GRCh38.p13_genomic.fna`


`pigeon filter classification.txt --isoforms sorted.gff`

Make Seurat compatible matrices for downstream analysis

`pigeon make-seurat --dedup dedup.fasta --group collapsed.group.txt -d seurat_out classification.filtered_lite_classification.txt`
