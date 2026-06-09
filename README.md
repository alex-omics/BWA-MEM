# BWA-MEM

A general-purpose WDL workflow for Illumina short-read alignment using BWA-MEM, with optional Picard duplicate marking and featureCounts quantification. Designed for Terra/Cromwell but compatible with any WDL 1.0 executor.
Suitable for: RNA-Seq, DNA-Seq, WGS, amplicon sequencing, ChIP-Seq, or any paired-end Illumina short-read data requiring alignment to a reference genome.

---

## Workflow structure

```
bwa_index          (once)   Build BWA index from reference FASTA
  └── bwa_mem      (scatter) Align PE FASTQs, sort, filter by MAPQ
        └── mark_duplicates  (scatter, optional) Picard MarkDuplicates
              └── featurecounts  (once, optional) featureCounts quantification
```

Alignment and duplicate marking scatter per sample. featureCounts runs once on all BAMs to produce a single merged count matrix. Both duplicate marking and featureCounts are optional and can be skipped via input parameters.

---

## Minimal usage (alignment only)

Supply only the required inputs — duplicate marking runs by default, featureCounts is skipped:

| Input | Type | Description |
|---|---|---|
| `read1_trimmed` | `Array[File]` | R1 FASTQ files (gzipped or plain) |
| `read2_trimmed` | `Array[File]` | R2 FASTQ files (gzipped or plain) |
| `sample_ids` | `Array[String]` | Sample identifiers (must match order of FASTQs) |
| `reference_fasta` | `File` | Reference genome FASTA |

---

## Full RNA-Seq usage

Additionally supply:

| Input | Type | Description |
|---|---|---|
| `reference_annotation` | `File?` | GFF3 or GTF annotation file. If provided, featureCounts runs. |

Key featureCounts parameters to set for RNA-Seq:

| Input | Default | Notes |
|---|---|---|
| `strandness` | `"2"` | `0`=unstranded, `1`=sense, `2`=antisense (reverse-stranded, e.g. Illumina Stranded Total RNA) |
| `feature_type` | `"CDS"` | `CDS` for bacteria, `exon` for eukaryotes |
| `attribute_type` | `"locus_tag"` | `locus_tag` for bacterial GFF, `gene_id` for eukaryotic GTF |
| `annotation_format` | `"GFF"` | `GFF` or `GTF` |

---

## All parameters

### Alignment

| Input | Type | Default | Description |
|---|---|---|---|
| `mapq_min` | `Int` | `20` | Minimum MAPQ to retain. `0`=no filter, `30`=strict |
| `mark_secondary` | `Boolean` | `true` | BWA `-M`: mark split hits as secondary (required for Picard compatibility) |
| `samtools_filter_flags` | `String` | `"0x904"` | `samtools view -F` flags. Excludes unmapped, secondary, supplementary. Set to `"0"` to disable filtering. |
| `bwa_extra_args` | `String` | `""` | Passthrough string for any additional BWA-MEM flags |

### Duplicate handling

| Input | Type | Default | Description |
|---|---|---|---|
| `skip_markdup` | `Boolean` | `false` | Skip Picard MarkDuplicates entirely (e.g. amplicon data) |
| `remove_duplicates` | `Boolean` | `false` | Remove duplicates rather than mark them. Not recommended for RNA-Seq. |

### featureCounts

| Input | Type | Default | Description |
|---|---|---|---|
| `strandness` | `String` | `"2"` | Library strandness |
| `feature_type` | `String` | `"CDS"` | GFF/GTF feature type to count |
| `attribute_type` | `String` | `"locus_tag"` | GFF/GTF attribute to use as gene ID |
| `annotation_format` | `String` | `"GFF"` | `GFF` or `GTF` |
| `paired_end_counting` | `Boolean` | `true` | Enable paired-end counting mode (`-p`) |
| `require_both_mates` | `Boolean` | `true` | Require both mates to map (`-B`) |
| `count_chimeric` | `Boolean` | `false` | Count chimeric read pairs (discouraged) |
| `ignore_duplicates` | `Boolean` | `true` | Ignore Picard-flagged duplicates (`--ignoreDup`) |
| `fraction_counting` | `Boolean` | `false` | Fractional counting for multi-mappers (`--fraction`) |
| `min_overlap` | `Int` | `1` | Minimum bases overlapping a feature (`--minOverlap`) |
| `output_prefix` | `String` | `"counts"` | Prefix for count matrix output filename |

### Compute resources

All resource parameters are independently overridable per task:

| Input | Default |
|---|---|
| `index_cpu` / `index_mem_gb` / `index_disk_gb` | `2` / `8` / `50` |
| `align_cpu` / `align_mem_gb` / `align_disk_gb` | `16` / `32` / `200` |
| `markdup_cpu` / `markdup_mem_gb` / `markdup_disk_gb` | `4` / `16` / `200` |
| `count_cpu` / `count_mem_gb` / `count_disk_gb` | `16` / `32` / `200` |

### Docker images

All images are overridable for version pinning or air-gapped environments:

| Input | Default |
|---|---|
| `bwa_docker` | `staphb/bwa:0.7.19` |
| `picard_docker` | `staphb/picard:3.1.0` |
| `subread_docker` | `biocontainers/subread:2.0.3--h9f5acd7_0` |

---

## Outputs

### Always produced (per sample)

| Output | Type | Description |
|---|---|---|
| `aligned_bams` | `Array[File]` | Sorted, filtered BAMs |
| `aligned_bais` | `Array[File]` | BAM indexes |
| `flagstats` | `Array[File]` | samtools flagstat per sample |
| `bwa_versions` | `Array[String]` | BWA version strings |
| `samtools_versions` | `Array[String]` | samtools version strings |

### If `skip_markdup=false` (default)

| Output | Type | Description |
|---|---|---|
| `markdup_bams` | `Array[File?]` | Duplicate-marked BAMs |
| `markdup_bais` | `Array[File?]` | BAM indexes |
| `dup_metrics` | `Array[File?]` | Picard duplication metrics per sample |

### If `reference_annotation` provided

| Output | Type | Description |
|---|---|---|
| `count_matrix` | `File?` | featureCounts merged count matrix (features × samples) |
| `featurecounts_summary` | `File?` | featureCounts assignment summary |

---

## Terra usage

1. Import this workflow into your Terra workspace via Dockstore or GitHub URL.
2. In the workflow configuration, select **"Run workflow(s) with inputs defined by data table"**.
3. Map `read1_trimmed` → `this.read1_trimmed`, `read2_trimmed` → `this.read2_trimmed`, `sample_ids` → `this.sample_id`.
4. Supply `reference_fasta` (and optionally `reference_annotation`) as workspace-level file paths.
5. Map desired outputs back to `this.aligned_bam`, `this.markdup_bam`, etc.

---

## Repository structure

```
illumina-bwamem-align-count/
├── .dockstore.yml
├── README.md
└── workflows/
    └── BWA_MEM_align.wdl
```

---

## Citation / acknowledgements

BWA-MEM: Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. arXiv:1303.3997.

Picard: http://broadinstitute.github.io/picard/

featureCounts: Liao Y, Smyth GK, Shi W. (2014) featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics 30(7):923-30.
