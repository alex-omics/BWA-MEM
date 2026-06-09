version 1.0

# =============================================================================
# illumina_bwamem_align_count.wdl
#
# General-purpose Illumina short-read alignment workflow using BWA-MEM.
# Optionally performs Picard duplicate marking and featureCounts quantification.
#
# Suitable for: RNA-Seq, DNA-Seq, WGS, amplicon, ChIP-Seq, or any paired-end
# Illumina short-read data requiring alignment to a reference genome.
#
# Minimal usage (alignment only):
#   Provide: read1_trimmed, read2_trimmed, sample_ids, reference_fasta
#   All other inputs have sensible defaults or are optional.
#
# Full RNA-Seq usage:
#   Additionally provide: reference_annotation
#   Set: strandness, feature_type, attribute_type, annotation_format
#
# Author: aa2084 (MGH)
# =============================================================================

workflow illumina_bwamem_align_count {

  meta {
    description: "General-purpose Illumina short-read alignment with BWA-MEM, optional Picard duplicate marking, and optional featureCounts quantification. Suitable for RNA-Seq, DNA-Seq, WGS, amplicon, and ChIP-Seq."
  }

  input {

    # -------------------------------------------------------------------------
    # Required: per-sample inputs (scatter)
    # Map to Terra data table columns e.g. this.read1_trimmed
    # -------------------------------------------------------------------------
    Array[File]   read1_trimmed
    Array[File]   read2_trimmed
    Array[String] sample_ids

    # -------------------------------------------------------------------------
    # Required: reference genome
    # -------------------------------------------------------------------------
    File   reference_fasta

    # -------------------------------------------------------------------------
    # Optional: annotation file for featureCounts quantification
    # If not provided, featureCounts is skipped and only BAMs are produced.
    # -------------------------------------------------------------------------
    File?  reference_annotation

    # -------------------------------------------------------------------------
    # Alignment parameters
    # -------------------------------------------------------------------------
    Int     mapq_min              = 20       # Min MAPQ to retain. 0=no filter, 20=default, 30=strict
    Boolean mark_secondary        = true     # BWA -M: mark split hits as secondary (Picard compat)
    String  samtools_filter_flags = "0x904"  # samtools view -F: excludes unmapped/secondary/supplementary
                                             # Set to "0" to disable all filtering
    String  bwa_extra_args        = ""       # Passthrough for any additional BWA-MEM flags

    # -------------------------------------------------------------------------
    # Duplicate handling
    # skip_markdup=true: skip MarkDuplicates entirely (e.g. amplicon data)
    # remove_duplicates=true: remove rather than mark (not recommended for RNA-Seq)
    # -------------------------------------------------------------------------
    Boolean skip_markdup          = false
    Boolean remove_duplicates     = false

    # -------------------------------------------------------------------------
    # featureCounts parameters (only used if reference_annotation is provided)
    # -------------------------------------------------------------------------
    String  strandness            = "2"        # 0=unstranded, 1=sense, 2=antisense (reverse-stranded)
    String  feature_type          = "CDS"      # CDS for bacteria, exon for eukaryotes
    String  attribute_type        = "locus_tag" # locus_tag (bacteria GFF), gene_id (eukaryote GTF)
    String  annotation_format     = "GFF"      # GFF or GTF
    Boolean paired_end_counting   = true       # -p: paired-end counting mode
    Boolean require_both_mates    = true       # -B: both mates must map (only if paired_end_counting)
    Boolean count_chimeric        = false      # if false, discard chimeric pairs (-C flag)
    Boolean ignore_duplicates     = true       # --ignoreDup: ignore Picard-flagged duplicates
    Boolean fraction_counting     = false      # --fraction: fractional counting for multi-mappers
    Int     min_overlap           = 1          # --minOverlap: min bases overlapping a feature
    String  output_prefix         = "counts"   # prefix for count matrix output filename

    # -------------------------------------------------------------------------
    # Compute resources (independently overridable per task)
    # -------------------------------------------------------------------------
    Int    index_cpu              = 2
    Int    index_mem_gb           = 8
    Int    index_disk_gb          = 50

    Int    align_cpu              = 16
    Int    align_mem_gb           = 32
    Int    align_disk_gb          = 200

    Int    markdup_cpu            = 4
    Int    markdup_mem_gb         = 16
    Int    markdup_disk_gb        = 200

    Int    count_cpu              = 16
    Int    count_mem_gb           = 32
    Int    count_disk_gb          = 200

    # -------------------------------------------------------------------------
    # Docker images (overridable for version pinning or air-gapped environments)
    # -------------------------------------------------------------------------
    String bwa_docker             = "staphb/bwa:0.7.19"
    String picard_docker          = "staphb/picard:3.1.0"
    String subread_docker         = "biocontainers/subread:2.0.3--h9f5acd7_0"
  }

  # ---------------------------------------------------------------------------
  # Task 1: Build BWA index (once per workflow run)
  # ---------------------------------------------------------------------------
  call bwa_index {
    input:
      reference_fasta = reference_fasta,
      cpu             = index_cpu,
      memory          = index_mem_gb,
      disk_size       = index_disk_gb,
      docker          = bwa_docker
  }

  # ---------------------------------------------------------------------------
  # Tasks 2-3: Align + optionally mark duplicates (scattered per sample)
  # ---------------------------------------------------------------------------
  scatter (i in range(length(sample_ids))) {

    call bwa_mem {
      input:
        read1                 = read1_trimmed[i],
        read2                 = read2_trimmed[i],
        samplename            = sample_ids[i],
        reference_fasta       = reference_fasta,
        index_amb             = bwa_index.index_amb,
        index_ann             = bwa_index.index_ann,
        index_bwt             = bwa_index.index_bwt,
        index_pac             = bwa_index.index_pac,
        index_sa              = bwa_index.index_sa,
        mapq_min              = mapq_min,
        mark_secondary        = mark_secondary,
        samtools_filter_flags = samtools_filter_flags,
        extra_args            = bwa_extra_args,
        cpu                   = align_cpu,
        memory                = align_mem_gb,
        disk_size             = align_disk_gb,
        docker                = bwa_docker
    }

    # Run MarkDuplicates only if skip_markdup=false (default)
    if (!skip_markdup) {
      call mark_duplicates {
        input:
          input_bam         = bwa_mem.bam,
          input_bai         = bwa_mem.bai,
          samplename        = sample_ids[i],
          remove_duplicates = remove_duplicates,
          cpu               = markdup_cpu,
          memory            = markdup_mem_gb,
          disk_size         = markdup_disk_gb,
          docker            = picard_docker
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Task 4: featureCounts — only runs if reference_annotation is provided
  # Uses markdup BAMs if available, otherwise falls back to aligned BAMs
  # ---------------------------------------------------------------------------
  if (defined(reference_annotation)) {

    # Select BAMs: prefer markdup output, fall back to aligned BAMs
    Array[File] bams_for_counting = if (!skip_markdup)
      then select_all(mark_duplicates.markdup_bam)
      else bwa_mem.bam

    call featurecounts {
      input:
        input_bams           = bams_for_counting,
        reference_annotation = select_first([reference_annotation]),
        strandness           = strandness,
        feature_type         = feature_type,
        attribute_type       = attribute_type,
        annotation_format    = annotation_format,
        paired_end           = paired_end_counting,
        require_both_mates   = require_both_mates,
        count_chimeric       = count_chimeric,
        ignore_duplicates    = ignore_duplicates,
        fraction_counting    = fraction_counting,
        min_overlap          = min_overlap,
        output_prefix        = output_prefix,
        cpu                  = count_cpu,
        memory               = count_mem_gb,
        disk_size            = count_disk_gb,
        docker               = subread_docker
    }
  }

  # ---------------------------------------------------------------------------
  # Workflow outputs
  # ---------------------------------------------------------------------------
  output {
    # Per-sample: always produced
    Array[File]   aligned_bams           = bwa_mem.bam
    Array[File]   aligned_bais           = bwa_mem.bai
    Array[File]   flagstats              = bwa_mem.flagstat
    Array[String] bwa_versions           = bwa_mem.bwa_version
    Array[String] samtools_versions      = bwa_mem.samtools_version

    # Per-sample: only if skip_markdup=false
    Array[File?]  markdup_bams           = mark_duplicates.markdup_bam
    Array[File?]  markdup_bais           = mark_duplicates.markdup_bai
    Array[File?]  dup_metrics            = mark_duplicates.metrics_file

    # Workflow-level: only if reference_annotation provided
    File?         count_matrix           = featurecounts.count_matrix
    File?         featurecounts_summary  = featurecounts.summary
  }
}


# =============================================================================
# TASK: bwa_index
# Builds BWA index from reference FASTA. Runs once per workflow execution.
# =============================================================================
task bwa_index {
  input {
    File   reference_fasta
    Int    cpu       = 2
    Int    memory    = 8
    Int    disk_size = 50
    String docker    = "staphb/bwa:0.7.19"
  }

  command <<<
    set -euo pipefail
    date | tee DATE
    echo "BWA $(bwa 2>&1 | grep Version)" | tee BWA_VERSION

    cp ~{reference_fasta} reference.fasta
    bwa index reference.fasta
  >>>

  output {
    File   index_amb   = "reference.fasta.amb"
    File   index_ann   = "reference.fasta.ann"
    File   index_bwt   = "reference.fasta.bwt"
    File   index_pac   = "reference.fasta.pac"
    File   index_sa    = "reference.fasta.sa"
    String bwa_version = read_string("BWA_VERSION")
  }

  runtime {
    docker:      docker
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB"
    preemptible: 1
    maxRetries:  2
  }
}


# =============================================================================
# TASK: bwa_mem
# Aligns paired-end FASTQs, sorts, and filters output BAM.
# =============================================================================
task bwa_mem {
  input {
    File    read1
    File    read2
    String  samplename
    File    reference_fasta
    File    index_amb
    File    index_ann
    File    index_bwt
    File    index_pac
    File    index_sa
    Int     mapq_min              = 20
    Boolean mark_secondary        = true
    String  samtools_filter_flags = "0x904"
    String  extra_args            = ""
    Int     cpu                   = 16
    Int     memory                = 32
    Int     disk_size             = 200
    String  docker                = "staphb/bwa:0.7.19"
  }

  command <<<
    set -euo pipefail
    date | tee DATE
    echo "BWA $(bwa 2>&1 | grep Version)" | tee BWA_VERSION
    samtools --version | head -n1 | tee SAMTOOLS_VERSION

    # Localize reference and all index files together
    cp ~{reference_fasta} reference.fasta
    cp ~{index_amb} reference.fasta.amb
    cp ~{index_ann} reference.fasta.ann
    cp ~{index_bwt} reference.fasta.bwt
    cp ~{index_pac} reference.fasta.pac
    cp ~{index_sa}  reference.fasta.sa

    # Input read count QC
    if [[ "~{read1}" == *".gz" ]]; then
      cat_cmd="zcat"
    else
      cat_cmd="cat"
    fi
    echo "R1 input reads: $("${cat_cmd}" "~{read1}" | awk 'NR%4==1' | wc -l)"
    echo "R2 input reads: $("${cat_cmd}" "~{read2}" | awk 'NR%4==1' | wc -l)"

    # Build optional flag arrays — bash arrays safely handle empty/absent flags
    # without the word-splitting risk of unquoted shell variables (SC2086)
    bwa_optional_flags=()
    if [ "~{mark_secondary}" == "true" ]; then
      bwa_optional_flags+=("-M")
    fi

    # extra_args passthrough — split on whitespace into array elements
    extra_args_array=()
    if [ -n "~{extra_args}" ]; then
      read -r -a extra_args_array <<< "~{extra_args}"
    fi

    samtools_filter_flags=()
    if [ "~{samtools_filter_flags}" != "0" ]; then
      samtools_filter_flags+=("-F" "~{samtools_filter_flags}")
    fi

    # Align, sort, filter
    bwa mem \
      -t ~{cpu} \
      "${bwa_optional_flags[@]+"${bwa_optional_flags[@]}"}" \
      -R "@RG\tID:~{samplename}\tSM:~{samplename}\tPL:ILLUMINA\tLB:~{samplename}\tPU:~{samplename}" \
      "${extra_args_array[@]+"${extra_args_array[@]}"}" \
      reference.fasta \
      "~{read1}" \
      "~{read2}" | \
    samtools sort \
      -@ ~{cpu} \
      -o ~{samplename}.sorted.bam -

    samtools view \
      -@ ~{cpu} \
      "${samtools_filter_flags[@]+"${samtools_filter_flags[@]}"}" \
      -q ~{mapq_min} \
      -b \
      -o ~{samplename}.sorted.filtered.bam \
      ~{samplename}.sorted.bam

    samtools index ~{samplename}.sorted.filtered.bam
    samtools flagstat ~{samplename}.sorted.filtered.bam | tee ~{samplename}.flagstat.txt
  >>>

  output {
    File   bam              = "~{samplename}.sorted.filtered.bam"
    File   bai              = "~{samplename}.sorted.filtered.bam.bai"
    File   flagstat         = "~{samplename}.flagstat.txt"
    String bwa_version      = read_string("BWA_VERSION")
    String samtools_version = read_string("SAMTOOLS_VERSION")
  }

  runtime {
    docker:      docker
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB"
    preemptible: 1
    maxRetries:  2
  }
}


# =============================================================================
# TASK: mark_duplicates
# Marks (or optionally removes) PCR/optical duplicates using Picard.
# =============================================================================
task mark_duplicates {
  input {
    File    input_bam
    File    input_bai
    String  samplename
    Boolean remove_duplicates = false
    Int     cpu               = 4
    Int     memory            = 16
    Int     disk_size         = 200
    String  docker            = "staphb/picard:3.1.0"
  }

  command <<<
    set -euo pipefail
    date | tee DATE
    picard MarkDuplicates --version 2>&1 | tee PICARD_VERSION || true

    picard MarkDuplicates \
      I=~{input_bam} \
      O=~{samplename}.markdup.bam \
      M=~{samplename}.dup_metrics.txt \
      REMOVE_DUPLICATES=~{remove_duplicates} \
      ASSUME_SORT_ORDER=coordinate \
      CREATE_INDEX=true \
      VALIDATION_STRINGENCY=LENIENT

    # Normalize index filename to .bam.bai convention
    if [ -f "~{samplename}.markdup.bai" ]; then
      mv ~{samplename}.markdup.bai ~{samplename}.markdup.bam.bai
    fi
  >>>

  output {
    File   markdup_bam  = "~{samplename}.markdup.bam"
    File   markdup_bai  = "~{samplename}.markdup.bam.bai"
    File   metrics_file = "~{samplename}.dup_metrics.txt"
  }

  runtime {
    docker:      docker
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB"
    preemptible: 1
    maxRetries:  2
  }
}


# =============================================================================
# TASK: featurecounts
# Counts reads per genomic feature across all samples simultaneously.
# Produces a single merged count matrix (features x samples).
# Only called when reference_annotation is provided.
# =============================================================================
task featurecounts {
  input {
    Array[File] input_bams
    File        reference_annotation
    String      strandness          = "2"
    String      feature_type        = "CDS"
    String      attribute_type      = "locus_tag"
    String      annotation_format   = "GFF"
    Boolean     paired_end          = true
    Boolean     require_both_mates  = true
    Boolean     count_chimeric      = false
    Boolean     ignore_duplicates   = true
    Boolean     fraction_counting   = false
    Int         min_overlap         = 1
    String      output_prefix       = "counts"
    Int         cpu                 = 16
    Int         memory              = 32
    Int         disk_size           = 200
    String      docker              = "biocontainers/subread:2.0.3--h9f5acd7_0"
  }

  command <<<
    set -euo pipefail
    date | tee DATE
    featureCounts -v 2>&1 | grep -i "featurecounts" | tee FEATURECOUNTS_VERSION || true

    # Build optional flag arrays to avoid SC2086 word-splitting issues
    fc_flags=()
    if [ "~{paired_end}" == "true" ]; then fc_flags+=("-p"); fi
    if [ "~{paired_end}" == "true" ] && [ "~{require_both_mates}" == "true" ]; then fc_flags+=("-B"); fi
    if [ "~{count_chimeric}" == "false" ]; then fc_flags+=("-C"); fi
    if [ "~{ignore_duplicates}" == "true" ]; then fc_flags+=("--ignoreDup"); fi
    if [ "~{fraction_counting}" == "true" ]; then fc_flags+=("--fraction"); fi

    featureCounts \
      -T ~{cpu} \
      "${fc_flags[@]+"${fc_flags[@]}"}" \
      -s ~{strandness} \
      -t ~{feature_type} \
      -g ~{attribute_type} \
      -F ~{annotation_format} \
      --minOverlap ~{min_overlap} \
      -a ~{reference_annotation} \
      -o ~{output_prefix}_matrix.txt \
      ~{sep=" " input_bams}
  >>>

  output {
    File count_matrix = "~{output_prefix}_matrix.txt"
    File summary      = "~{output_prefix}_matrix.txt.summary"
  }

  runtime {
    docker:      docker
    memory:      memory + " GB"
    cpu:         cpu
    disks:       "local-disk " + disk_size + " SSD"
    disk:        disk_size + " GB"
    preemptible: 0
    maxRetries:  2
  }
}
