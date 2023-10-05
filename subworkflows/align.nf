process call_paftools {
    // Covert gtf/gff3 files to bed12

    label "singlecell"
    cpus 1
    input:
        path "ref_genes.gtf"
    output:
        path "ref_genes.bed", emit: ref_genes_bed
    """
    paftools.js gff2bed -j ref_genes.gtf > ref_genes.bed
    """
}

process get_chrom_sizes{
    // cut -f1,2 extracts the 1st and 2nd columns of the genome index file (.fai)
    // sort -V sorts using version sort (-V) that is aware of numbers within the text (i.e. chromosome names)
    label "singlecell"
    cpus 1
    input:
        path "ref_genome.fai"
    output:
        path 'chr_sizes', emit: ref_chrom_sizes
    """
    cut -f1,2 ref_genome.fai | sort -V > chr_sizes
    """
}

process align_to_ref {
    // Map long-read sequences to reference genome using minimap2
    // -ax splice is used to aling RNA-seq data that contain spliced elements
    // -uf specifies that input file must be in FASTQ format and output in SAM format
    // --secondary=no indicates that no secondary alignments should be generated
    // --MD option used to include the MD tags (mismatches and deletions) in SAM output
    // -t specifies the number of CPU cores to use for the alignment
    // --junc-bed provides a BED file to guide the spliced alignments process
    // $params.resources_mm2_flags variable includes additional flags or parameters that are passed to minimap2
    // samtools view -b convert the SAM files into BAM format
    // --no-PG tag indicates not to add the @PG line to the header of the BAM file
    // -t ref_chrom_sizes specifies a reference chromosome sizes file that is used to set the header information for the BAM file
    // - input comes from standard input (piped fro previous commands)
    // samtools sort sorts the BAM file by genomic coordinates
    // -@ specifies the number of CPU threads should be used for sorting
    // samtools index creates an index file for the sorted BAM file
    
    label "singlecell"
    cpus params.resources_mm2_max_threads
    input:
        tuple val(sample_id),
              path("reads.fastq")
        path "ref_genome.fasta"
        path "ref_genes.bed"
        path "ref_chrom_sizes.tsv"
    output:
        tuple val(sample_id), 
            path("*sorted.bam"), 
            path("*sorted.bam.bai"), 
            emit: bam_sort
    """
     minimap2 -ax splice -uf --secondary=no --MD -t $task.cpus \
      --junc-bed ref_genes.bed $params.resources_mm2_flags  \
      ref_genome.fasta reads.fastq* \
        | samtools view -b --no-PG -t ref_chrom_sizes - \
        | samtools sort -@ 2 --no-PG  - > "${sample_id}_sorted.bam"
    samtools index -@ ${task.cpus} "${sample_id}_sorted.bam"
    """
}


// workflow module
workflow align {
    take:
        stranded_fq
        ref_genome
        ref_genome_idx
        ref_genes_gtf
    main:
        call_paftools(ref_genes_gtf)
        get_chrom_sizes(ref_genome_idx)
        align_to_ref(
            stranded_fq.groupTuple(),
            ref_genome,
            call_paftools.out.ref_genes_bed,
            get_chrom_sizes.out.ref_chrom_sizes)
    emit:
        bam_sort = align_to_ref.out.bam_sort
}
