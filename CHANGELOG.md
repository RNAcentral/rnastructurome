# nf-core/rnastructurome: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [unreleased]

**Behaviour change:** treated groups whose `sample_group` shares no leading token with any control now get no control; the reference-wide lone-control fallback is gone.

### `Added`

- `--cutadapt_nextseq_trim` passes Cutadapt `--nextseq-trim` with the given 3′ quality cutoff, stripping the high-quality poly-G tails that 2-colour (NextSeq/NovaSeq) chemistry leaves behind and quality trimming misses.
- `--cutadapt_discard_untrimmed` drops reads in which no adapter was found, for short-insert libraries (e.g. tRNA, miRNA) where every genuine read runs into the adapter. It cannot be combined with `--cutadapt_quality_only`.
- `--rfcount_map_max_coverage` passes rf-count's `--max-coverage` for MaP samples, downsampling each transcript above the cap (minimum `1000`), so a transcript at extreme depth (e.g. a viral genome at ~10<sup>6</sup>×) counts in hours rather than days. Off by default. Ignored with `--count_genome true`.

### `Changed`

- Fuzzy control pairing (`--fuzzy_untreated_pairing`) now picks the untreated or denatured control whose `sample_group` shares the longest leading token prefix with the treated group, rather than matching on the first token only, so each arm of a design such as `HFF_uninfected` / `HFF_HCMV_05hpi` / `HFF_HCMV_72hpi` gets its own control instead of aborting as ambiguous. A lone control is reused across an arm's replicates, and denatured controls are resolved the same way.
- `rf-count` now counts only the references that have alignments, instead of running one `samtools view` per FASTA entry; whole-transcriptome runs (~250k human transcripts) previously never finished. RC files therefore list only covered references. The `codon` profile now gives the genome-route `rf-count` task the same resources as the transcriptome-route one.

### `Fixed`

- Fuzzy control pairing only borrows a control from the same reference, so samples from different organisms can no longer pair.
- A `sample_group`/`replicate` group whose samples come from different organisms now stops with an error instead of pairing them.
- A `sample_group`/`replicate` group with more than one untreated or denatured sample now stops with a clear error instead of crashing; give runs of the same control the same sample name so they are concatenated.
- `rf-count` no longer dies with "Unable to extract" when any of the first BAM records lacks an MD tag (an unmapped mate is enough): it now receives a mapped-only, MD-tagged, indexed BAM.
- `rf-count` no longer hangs at 100% CPU on references containing IUPAC ambiguity codes (e.g. `Y`, `R`); these are masked to `N` before counting, so those positions report no reactivity.
- R2DT no longer fails on a missing `.command.env` when no sequences are extracted for drawing.
- `--rfnorm_norm_method 1` (2-8% normalisation) is accepted again; the launch guard and schema rejected it while `docs/usage.md` documented it.
- The rf-fold flag letters quoted in the `nextflow_schema.json` descriptions of `rffold_unconstrained`, `rffold_vienna_no_lonely_pairs`, `rffold_vienna_constrained`, `rffold_vienna_max_bp_span`, `rffold_fold_constraint_file` and `rffold_dotplot` now match what the pipeline passes (`-i`, `-nlp`, `-hc`, `-md`, `-c`, `-dp`).
- `.bp` arc files list base pairs sorted by position instead of in rf-fold's arbitrary dot-plot order, so reruns give identical files.
- `count/rfcount_summary_all_samples.tsv` lists samples in name order instead of the order their `rf-count` tasks finished, so reruns give identical files.
- `rf-normfactor` no longer silently drops RC files whose sample name starts with a digit (e.g. `125ng_r1.rc`): `rf-rctools index` numifies such a bare filename into an argument index, so the file is now passed with a directory prefix.

### `Removed`

- `--rffold_vienna_bp_span` and `--rffold_unpaired_constraint_file`, which were declared but never used; rf-fold has no corresponding flags.

## v1.0.0 - [8 September 2026]

First release of nf-core/rnastructurome, which analyses chemical high-throughput RNA structure-probing data and predicts RNA secondary structures from it. The pipeline covers **SHAPE** and **DMS** chemistries read out by either the **RT-stop** or **mutational profiling (MaP)** principle, and runs in two modes depending on the reference: a **genome** route aligned with STAR and a **transcriptome** route aligned with Bowtie/Bowtie2.

### `Added`

- Single- and paired-end FASTQ input, with re-sequenced samples concatenated automatically.
- Read quality control with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) before and after adapter and quality trimming with [Cutadapt](https://cutadapt.readthedocs.io/).
- Optional UMI extraction and UMI-aware deduplication with [UMI-tools](https://umi-tools.readthedocs.io/); position-based duplicate removal with [SAMtools markdup](https://www.htslib.org/) is available but off by default.
- Reference genomes and annotations downloaded automatically from [Ensembl](https://www.ensembl.org/), falling back to [NCBI](https://www.ncbi.nlm.nih.gov/datasets/) for bacteria, viruses and other organisms Ensembl does not cover, or supplied directly as local files.
- Genome-route alignment with [STAR](https://github.com/alexdobin/STAR), and transcriptome-route alignment with [Bowtie](https://bowtie-bio.sourceforge.net/) for RT-stop and [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/) for MaP.
- BAM sorting, indexing and alignment statistics with [SAMtools](https://www.htslib.org/).
- Per-base reactivity counting with [rf-count](https://rnaframework-docs.readthedocs.io/en/latest/rf-count/), tallying mutations for MaP and RT-stops for RT-stop on transcript coordinates.
- An alternative genome-coordinate counting path with [rf-count-genome](https://rnaframework-docs.readthedocs.io/en/latest/rf-count-genome/), resolving library strandedness with [BEDOPS](https://bedops.readthedocs.io/) and [RSeQC](https://rseqc.sourceforge.net/) and extracting per-transcript reactivity with [rf-rctools](https://rnaframework-docs.readthedocs.io/en/latest/rf-rctools/).
- Reactivity normalisation with [rf-norm](https://rnaframework-docs.readthedocs.io/en/latest/rf-norm/), pairing treated, untreated and denatured samples automatically and selecting scoring and normalisation methods from the controls present.
- Replicate reproducibility QC with [rf-correlate](https://rnaframework-docs.readthedocs.io/en/latest/rf-correlate/), reporting pairwise Pearson and Spearman correlation of reactivity profiles.
- RNA secondary structure prediction across grouped replicates with [rf-fold](https://rnaframework-docs.readthedocs.io/en/latest/rf-fold/), using chemistry- and reagent-aware folding defaults.
- Reactivity, Shannon entropy and base-pair arc tracks in transcript and genome coordinates, generated with [rf-wiggle](https://rnaframework-docs.readthedocs.io/en/latest/rf-wiggle/) and [UCSC wigToBigWig](https://genome.ucsc.edu/).
- 2D structure diagrams coloured by reactivity, drawn with [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) RNAplot and, where a template model exists for the RNA type, with [R2DT](https://github.com/RNAcentral/R2DT) alongside them. R2DT is container-only, so `-profile conda` draws with ViennaRNA alone.
- RMDB-compatible RDAT export combining per-transcript reactivity and structure.
- Optional folding calibration against known structures with [rf-jackknife](https://rnaframework-docs.readthedocs.io/en/latest/rf-jackknife/), which tunes the slope and intercept passed to rf-fold.
- Optional structural-element extraction of high-confidence, low-reactivity and low-Shannon motifs with [rf-structextract](https://rnaframework-docs.readthedocs.io/en/latest/rf-structextract/).
- Optional structure-accuracy evaluation with [rf-eval](https://rnaframework-docs.readthedocs.io/en/latest/rf-eval/), reporting AUROC, DSCI and the unpaired coefficient against an automatic rotation baseline.
- Aggregated quality-control report with [MultiQC](http://multiqc.info/), including alignment, reactivity and replicate-correlation summaries.
- Test profiles for the genome route, the transcriptome route, a prokaryote run with jackknife calibration, and a full-size rice DMS-MaPseq dataset.

### `Fixed`

### `Dependencies`

### `Deprecated`
