# HiFi assembly workflow: contributor guidance

## Scope and layout

This is a Linux x86_64 Snakemake workflow for PacBio HiFi assembly with hifiasm
and/or LJA, OATK mitochondrial recovery, overlap-based duplicate purging,
optional RagTag scaffolding, and assembly QC. Illumina reads supply additional
mapping and k-mer evidence. MGA is an optional, independent first-phase branch.

- `environment.yaml`: the Snakemake/Python runtime, including pandas.
- `config/config.yaml`, `config/samples.tsv`: generic public templates.
- `workflow/Snakefile`: input validation, sample selections, includes, targets.
- `workflow/accessors.smk`: shared paths and assembly-stage selection.
- `workflow/rules/`, `workflow/envs/`, `workflow/scripts/`: execution and helpers.
- `docs/figures/`: editable release figures, exports, and drawing source.
- `tests/`: synthetic fixtures, DAG tests, and independent output validation.

Only included and targeted rules execute. MultiQC, TIDK, heterozygosity, aggregate
QC, and MGA are disabled by default. MGA ends at native polished FASTA and graphs;
it does not enter the core purge/RagTag/final/QC branches. PloidyFrost, Smudgeplot,
ploidy plotting, and the separate RagTag Chromeister include are experimental,
inactive modules. Do not enable them incidentally.
The obsolete `filter_hifiasm`/`filter_lja` rules and their unused environment
stubs have been retired. Keep mitochondrial integration in the supported final
assembly rules; do not restore the old 50% PAF-coverage shortcut.

## Configuration and data boundaries

1. Inspect Git status and applicable instructions before edits; preserve unrelated
   user changes. Existing projects may have their own `AGENTS.md` instructions.
2. Copy public templates to ignored `.local` files. Keep host-specific sample
   tables, reads, references, databases, environments, and results out of Git.
3. Use absolute paths and unique, safe sample IDs. The Snakefile creates and
   switches to `project_dir` during parsing, including dry-runs. Use a new project
   directory for validation; a CLI working-directory override does not cancel
   the workflow's `workdir` directive.
4. At least one of hifiasm/LJA is required per sample. Sample-sheet `run_mga`
   overrides the YAML fallback. Required read/reference paths are checked before
   DAG construction; `Ref` remains required when RagTag is disabled.
5. Code fixes do not repair historical results. Keep production data unchanged
   during reviews; requested reruns should use separate output locations unless
   replacement was explicitly requested. Record inputs, commands, versions,
   parameters, checksums, and failed attempts.
6. Distinguish raw contigs, purged primary/hap sequences, reference-guided
   scaffolds, mitochondrial integration, and later curation. Candidate final
   FASTAs are not automatically validated curated genomes. Transfer QC between
   assemblies only after checking sequence identity, not sample names alone.

## Implementation and scientific invariants

- Preserve sample, assembler, and processing-stage identity throughout the DAG.
  Disabling an assembler must not schedule it through downstream dependencies.
- Declare real artifacts, helper scripts, reference/database inputs, threads,
  and meaningful memory requirements. Give each concurrent job its own mutable
  workspace; separate diagnostic stderr from FASTA/GFF/report stdout.
- Completion markers follow successful content checks. Never manufacture empty
  figures, swallow tool failures, or convert untestable results into zero or a
  default biological/model assignment. Preserve failure diagnostics.
- Purging currently uses overlaps without read-depth classification. Keep hap
  outputs and review large removals, especially for hybrid/polyploid objectives.
  Overlap labels alone do not establish that removed sequence is dispensable.
- Mitochondrial replacement uses assembly-target PAF IDs and unions split/circular
  alignments. Preserve nuclear flanks and unique identifiers for appended records.
- Merqury uses separate read/assembly database names and an independent workspace
  for every sample/stage/assembler. Validate report labels and real plot outputs.
- GenomeScope2 currently pools HiFi and both Illumina mates. Extract actual model
  fit summaries; the model ranking does not establish cellular ploidy.
- nQuire compares fixed models using one preferred final BAM per sample. Empty,
  non-informative, malformed, or tied evidence must not become an arbitrary model.
- Heterozygosity is a depth-callable diploid-model estimate using matching final
  FASTA/BAMs, originally biallelic SNPs, consistent depth/MQ/BQ/allele filters,
  and explicit exclusions. It does not confirm biological heterozygosity/ploidy.
- TIDK reports exact-repeat candidates and end coordinates, not certified T2T
  assemblies or chromosome counts. Keep missingness and limitations explicit.
- Aggregate QC retains BUSCO lineage/marker context, infinite QV, reference/input
  hashes, paired-read denominators, and optional-module statuses. MultiQC is a
  separate supported-native-report route; exclude BUSCO run-directory aliases.

## Validation and release checks

Create the documented runtime with `conda env create --file environment.yaml`.
Use the validated runtime and per-rule environments. Environment provisioning and
source retrieval require network access; current Kaleido plotting requires
Chrome/Chromium. GenomeScope2's required R dependencies come from Conda; install
its pinned local source without fetching optional CRAN packages.

Run relevant regressions with `python -m unittest discover -s tests -v`.
Native-tool fixtures require the relevant executables; report skips honestly.
Exercise both core assemblers, each single-assembler selection, RagTag on/off,
and relevant optional flags for DAG changes. Placeholder inputs prove DAG
construction, not tool execution or scientific correctness. Use synthetic native
fixtures for caller/plotter behavior and isolated full-read runs when requested.

Validate YAML/TSV types, helper syntax, fresh environment resolution, meaningful
artifacts, and output provenance separately. A Conda solve does not validate pip
installation, source builds, external browsers, or biological output quality.
Keep release audit evidence in ignored `test_runs/`; it may contain local paths
and dataset-specific results and should not be published by default.

Regenerate figures using `python docs/figures/render_workflows.py`. Confirm their
connections still match active rules and optional-module scope. Before release,
check Git-visible files, ignore behavior, local links, generated figures, and
required untracked source files. Publishing/tagging is a separate user action.
