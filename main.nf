#!/usr/bin/env nextflow
/*
 * brainvar-star-2pass
 * -------------------
 * Faithful Nextflow re-implementation of
 *   bmrn-wwrd/CNS_Genomics_Pipelines :: code/snakemake-pipeline/rna-seq-star
 * (STAR 2.7.10b two-pass, bbduk trimming) as used for the BrainVar cohort.
 *
 * Deliverables are per-run SJ.out.tab + Log.final.out. BAMs are produced by
 * pass 2 (as in the reference pipeline) but are NOT published.
 *
 * Parameter provenance: every STAR / bbduk flag below is copied from the
 * reference Snakefile rules trim / starindex / star1 / sjmerge / star2.
 */

nextflow.enable.dsl = 2

params.cohort        = null          // ngn2 | mouse | rhesus
params.samplesheet   = null          // CSV: sample,run_accession,fastq_1,fastq_2,...
params.genome_url    = null          // bgzipped primary-assembly FASTA
params.gtf_url       = null          // bgzipped GTF
params.outdir        = null
params.minsj         = 3             // config.yaml params.minsj
params.imax          = 1000000       // config.yaml params.imax
params.glen          = 0             // config.yaml params.glen
params.sjover        = 6             // config.yaml params.sjover
params.sjdbover      = 100           // config.yaml params.sjdbover
params.kmer          = 23            // config.yaml params.kmer
params.mink          = 11            // config.yaml params.mink
params.hdist         = 1             // config.yaml params.hdist
params.keep_contigs  = null          // regex of contigs to KEEP (rhesus/mouse); null => repo grep -v
params.transcriptome = false         // repo emits TranscriptomeSAM for RSEM; we do not run RSEM

if( !params.cohort )      error "must supply --cohort"
if( !params.samplesheet ) error "must supply --samplesheet"
if( !params.outdir )      error "must supply --outdir"


/* ------------------------------------------------------------------ *
 *  References
 * ------------------------------------------------------------------ */

process FETCH_GENOME {
    tag "${params.cohort}"
    label 'net'
    storeDir "${params.outdir}/_refs"

    output:
    path "genome.fa.gz"

    script:
    """
    set -euo pipefail
    for a in 1 2 3 4 5 6 7 8; do
      if wget --continue --timeout=60 --tries=3 --no-verbose -O genome.fa.gz '${params.genome_url}' \\
         && gzip -t genome.fa.gz 2>/dev/null; then
        exit 0
      fi
      rm -f genome.fa.gz
      sleep \$(( a * 45 + RANDOM % 60 ))
    done
    echo "FATAL: genome download failed" >&2; exit 1
    """
}

process FETCH_GTF {
    tag "${params.cohort}"
    label 'net'
    storeDir "${params.outdir}/_refs"

    output:
    path "annotation.gtf.gz"

    script:
    """
    set -euo pipefail
    for a in 1 2 3 4 5 6 7 8; do
      if wget --continue --timeout=60 --tries=3 --no-verbose -O annotation.gtf.gz '${params.gtf_url}' \\
         && gzip -t annotation.gtf.gz 2>/dev/null; then
        exit 0
      fi
      rm -f annotation.gtf.gz
      sleep \$(( a * 45 + RANDOM % 60 ))
    done
    echo "FATAL: gtf download failed" >&2; exit 1
    """
}

/*
 * rule starindex — NOTE: the reference Snakefile builds the index WITHOUT
 * --sjdbGTFfile / --sjdbOverhang. Annotation enters at pass 2 via on-the-fly
 * sjdb insertion. We reproduce that exactly.
 */
process STAR_INDEX {
    tag "${params.cohort}"
    label 'star_big'

    input:
    path fasta

    output:
    path "starindex"

    script:
    """
    gzip -t ${fasta}
    zcat ${fasta} > genome.fa
    mkdir -p starindex
    STAR --runThreadN ${task.cpus} \\
      --runMode genomeGenerate --genomeDir starindex \\
      --genomeFastaFiles genome.fa \\
      --outTmpDir \$PWD/_startmp
    rm -f genome.fa
    """
}


/* ------------------------------------------------------------------ *
 *  Per-run: fetch -> trim -> pass 1
 * ------------------------------------------------------------------ */

process FETCH_FASTQ {
    tag "$run"
    label 'net'
    maxForks params.fetch_forks

    input:
    tuple val(run), val(url1), val(url2)

    output:
    tuple val(run), path("${run}_R1.fastq.gz"), path("${run}_R2.fastq.gz")

    script:
    // ENA ftp:// paths are served identically over https.
    // EBI refuses connections (not 4xx — TCP refused) when too many streams
    // originate from one NAT IP, so: download the mates SERIALLY, with long
    // randomised backoff, and cap concurrency via params.fetch_forks.
    def u1 = url1.replaceFirst(/^ftp:\/\//, 'https://')
    def u2 = url2.replaceFirst(/^ftp:\/\//, 'https://')
    """
    set -euo pipefail

    fetch () {
      local url="\$1" out="\$2"
      for attempt in 1 2 3 4 5 6 7 8 9 10; do
        if wget --continue --timeout=60 --waitretry=30 --tries=3 \
                --no-verbose -O "\$out" "\$url"; then
          if gzip -t "\$out" 2>/dev/null; then
            echo "OK \$out on attempt \$attempt" >&2
            return 0
          fi
          echo "gzip check failed for \$out (attempt \$attempt), refetching" >&2
          rm -f "\$out"
        fi
        # exponential-ish backoff with jitter, capped
        local wait=\$(( attempt * 45 + RANDOM % 60 ))
        echo "attempt \$attempt failed for \$url; sleeping \${wait}s" >&2
        sleep "\$wait"
      done
      echo "FATAL: could not download \$url after 10 attempts" >&2
      return 1
    }

    # stagger task starts so a burst of tasks does not hit EBI simultaneously
    sleep \$(( RANDOM % 45 ))

    fetch '${u1}' ${run}_R1.fastq.gz
    sleep \$(( 5 + RANDOM % 15 ))
    fetch '${u2}' ${run}_R2.fastq.gz

    test -s ${run}_R1.fastq.gz
    test -s ${run}_R2.fastq.gz
    """
}

/* rule trim */
process BBDUK_TRIM {
    tag "$run"
    label 'trim'
    publishDir "${params.outdir}/qc/bbduk", mode: 'copy', pattern: "*_bbduk_stat.txt"

    input:
    tuple val(run), path(fq1), path(fq2)

    output:
    tuple val(run), path("t_${run}_R1.fastq.gz"), path("t_${run}_R2.fastq.gz"), emit: reads
    path "${run}_bbduk_stat.txt", emit: stat

    script:
    """
    ADAPTERS=\$(find / -name adapters.fa -path '*bbmap*' 2>/dev/null | head -1)
    if [ -z "\$ADAPTERS" ]; then echo "bbmap adapters.fa not found" >&2; exit 1; fi
    echo "using adapters: \$ADAPTERS" >&2

    bbduk.sh -Xmx${task.memory.toGiga() - 2}g in1=${fq1} in2=${fq2} \\
      out1=t_${run}_R1.fastq.gz out2=t_${run}_R2.fastq.gz t=${task.cpus} \\
      ref=\$ADAPTERS stats=${run}_bbduk_stat.txt ktrim=r k=${params.kmer} \\
      mink=${params.mink} hdist=${params.hdist} tpe tbo
    """
}

/* rule star1 — splice junctions only, no BAM */
process STAR_PASS1 {
    tag "$run"
    label 'star'
    publishDir "${params.outdir}/pass1_sj", mode: 'copy', pattern: "*_SJ.out.tab"

    input:
    tuple val(run), path(fq1), path(fq2)
    path index

    output:
    path "${run}_SJ.out.tab", emit: sj
    path "${run}_Log.final.out", emit: log

    script:
    """
    STAR --runThreadN ${task.cpus} --genomeDir ${index} \\
      --readFilesIn ${fq1} ${fq2} \\
      --outFileNamePrefix ${run}_ \\
      --outSAMtype None --alignIntronMax ${params.imax} \\
      --scoreGenomicLengthLog2scale ${params.glen} \\
      --alignSJoverhangMin ${params.sjover} \\
      --readFilesCommand zcat \\
      --outTmpDir \$PWD/_startmp
    test -s ${run}_SJ.out.tab
    """
}


/* ------------------------------------------------------------------ *
 *  Cohort-level junction database
 * ------------------------------------------------------------------ */

/*
 * rule sjmerge. Reference implementation:
 *   cat *_SJ.out.tab | grep -v -E '^G|J|K|M|chrM|chrU|HLA' | sort | uniq -c \
 *     | awk '{OFS="\t"}{if ($1 >= minsj) {print $2,...,$10}}'
 *
 * The repo regex is GRCh38/UCSC-specific. Per species we substitute an
 * equivalent scaffold filter (params.keep_contigs) but keep the
 * count/threshold logic byte-identical.
 */
process SJ_MERGE {
    tag "${params.cohort}"
    label 'small'
    publishDir "${params.outdir}/sjdb", mode: 'copy'

    input:
    path sjtabs

    output:
    path "splice-junctions-filtered.tsv", emit: sjdb
    path "sjmerge_counts.txt",            emit: counts

    script:
    def filt = params.keep_contigs
        ? "awk -F'\\t' '\$1 ~ /^(${params.keep_contigs})\$/'"
        : "grep -v -E '^G|J|K|M|chrM|chrU|HLA'"
    """
    cat ${sjtabs} | ${filt} | sort | uniq -c \\
      | awk '{OFS="\\t"}{if (\$1 >= ${params.minsj}) {print \$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}}' \\
      > splice-junctions-filtered.tsv

    # diagnostics: total observed vs retained
    TOTAL=\$(cat ${sjtabs} | ${filt} | cut -f1-4 | sort -u | wc -l)
    KEPT=\$(cut -f1-4 splice-junctions-filtered.tsv | sort -u | wc -l)
    LINES=\$(wc -l < splice-junctions-filtered.tsv)
    {
      echo "cohort\t${params.cohort}"
      echo "n_input_sj_files\t\$(echo ${sjtabs} | wc -w)"
      echo "unique_junctions_after_contig_filter\t\$TOTAL"
      echo "rows_in_filtered_sjdb\t\$LINES"
      echo "unique_junctions_in_filtered_sjdb\t\$KEPT"
      echo "minsj\t${params.minsj}"
    } > sjmerge_counts.txt
    test -s splice-junctions-filtered.tsv
    """
}


/* ------------------------------------------------------------------ *
 *  Pass 2
 * ------------------------------------------------------------------ */

/* rule star2 */
process STAR_PASS2 {
    tag "$run"
    label 'star_big'
    publishDir "${params.outdir}/pass2_sj",  mode: 'copy', pattern: "*_SJ.out.tab"
    publishDir "${params.outdir}/logs",      mode: 'copy', pattern: "*_Log.final.out"
    publishDir "${params.outdir}/genecounts", mode: 'copy', pattern: "*_ReadsPerGene.out.tab"

    input:
    tuple val(run), path(fq1), path(fq2)
    path index
    path sjdb
    path gtf

    output:
    path "${run}_SJ.out.tab",           emit: sj
    path "${run}_Log.final.out",        emit: log
    path "${run}_ReadsPerGene.out.tab", emit: counts

    script:
    def quant = params.transcriptome ? 'TranscriptomeSAM GeneCounts' : 'GeneCounts'
    """
    zcat ${gtf} > annotation.gtf

    STAR --runThreadN ${task.cpus} --genomeDir ${index} \\
      --readFilesIn ${fq1} ${fq2} \\
      --outFileNamePrefix ${run}_ \\
      --sjdbFileChrStartEnd ${sjdb} --sjdbGTFfile annotation.gtf \\
      --sjdbOverhang ${params.sjdbover} \\
      --outSAMtype BAM SortedByCoordinate --outFilterType BySJout \\
      --alignIntronMax ${params.imax} \\
      --scoreGenomicLengthLog2scale ${params.glen} \\
      --alignSJoverhangMin ${params.sjover} \\
      --quantMode ${quant} \\
      --outSAMunmapped Within --outSAMattributes All \\
      --readFilesCommand zcat \\
      --outTmpDir \$PWD/_startmp

    test -s ${run}_SJ.out.tab
    rm -f *.bam annotation.gtf
    """
}


/* ------------------------------------------------------------------ *
 *  Workflow
 * ------------------------------------------------------------------ */

workflow {

    rows = Channel
        .fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { r -> tuple(r.run_accession, r.fastq_1, r.fastq_2) }

    genome = FETCH_GENOME()
    gtf    = FETCH_GTF()
    index  = STAR_INDEX(genome)

    raw    = FETCH_FASTQ(rows)
    trim   = BBDUK_TRIM(raw)

    p1     = STAR_PASS1(trim.reads, index)

    sjdb   = SJ_MERGE(p1.sj.collect())

    STAR_PASS2(trim.reads, index, sjdb.sjdb, gtf)
}

workflow.onComplete {
    log.info "cohort=${params.cohort} status=${workflow.success ? 'OK' : 'FAILED'} duration=${workflow.duration}"
}
