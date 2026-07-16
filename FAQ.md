# FAQ — Lessons & Limitations from Real Analyses

[← Back to README](README.md)

Non-obvious issues and design decisions distilled from production Visium HD runs
(the intestinal Duodenum/Ileum project). Each answer includes a concrete example.

## Environment & HPC (Minerva / LSF)

**Q: My R packages won't load on the cluster — I get ABI / symbol errors. What's wrong?**
Always `module purge && module load R/4.4.1` and set
`.libPaths("/sc/arion/projects/naiklab/ikjot/R")`. That library is compiled
against the Rocky9 R/4.4.1 module; **conda R or any other R build causes ABI
mismatches and package load failures.** *Example:* switching to a conda R env to
"save time" broke `Seurat`/`spacexr` loading with cryptic symbol errors — reverting
to the module fixed it instantly. Never use conda R for project scripts.

**Q: My single-core job (e.g. RCTD) segfaults with no R error message. Why?**
OpenBLAS spawns one thread **per physical core** by default. When a step requests a
single core (`create.RCTD(max_cores = 1)`), the extra threads exhaust `RLIMIT_NPROC`
and the process segfaults silently. **Fix:** pin all BLAS/OMP backends to one thread
(`export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLAS_NUM_THREADS=1`).
*Example:* RCTD on a 372k-spot 8µm sample died with no traceback until these caps
were added — see `submit_template.sh`. This is the #1 "job died mysteriously" cause.

**Q: How much memory / walltime do I actually need?**
Profiles that worked in practice: **8µm RCTD** (~372k spots) needed **300 GB/core × 2
(600 GB)** after OOM-ing at 200 GB; **16µm RCTD** ran at ~200 GB/core; both on the
`premium` queue with **144 h** walltime and `span[hosts=1]`. `rusage[mem=...]` is **per
core**, not total. Start high on the largest sample and trim down.

**Q: `.Rmd` or `.R` on the cluster?**
Run heavy compute as a plain **`.R` via `bsub`** (`Rscript step.R`); render `.Rmd`
only for lighter reporting. *Caution:* keep the two in sync — a stale `.R` diverging
from its notebook silently ran old logic in a past project.

## RCTD Deconvolution

**Q: Can I use one shared reference for all samples/tissues?**
Prefer **tissue/group-specific references**, and understand the trade-off — this is a
real **limitation**, not just a setting. If a cell type is structurally absent from a
tissue but left in that tissue's reference, RCTD can assign it spuriously; but if you
*exclude* it, any cross-tissue comparison of that type becomes **"present in reference"
vs "structurally forced to zero"** — an artifact, not biology. *Example:*
`Mid_villus_proximal_enterocytes` is excluded from the Ileum reference and
`Mid_villus_distal_enterocytes` from the Duodenum reference; we therefore **drop those
types entirely** from cross-tissue composition plots and tests rather than compare them.

**Q: RCTD errors on my labels or counts.**
Three recurring fixes: (a) RCTD **prohibits `/` in cell-type labels** → sanitize to `-`;
(b) merged **SCT objects store non-integer RNA counts** → `round()` and coerce to
`dgCMatrix` before building the reference; (c) drop rare types (**min ~25 cells**) and
special labels (`Unknown`, `Mix cluster remove`). *Example:* a `TA_cells / Early_progenitors`
label crashed `Reference()` until the slash was replaced.

**Q: Do I have to recompute RCTD every run?**
No — **cache per-sample RCTD objects to RDS and reload**, recomputing only when the file
is missing. A full re-run of a 4-sample 8µm deconvolution is many hours; caching makes
iterating on the downstream plots cheap.

## Reference & Annotation Hygiene

**Q: Cell types duplicate, split, or lose their color in plots. Why?**
Label whitespace/Unicode inconsistencies in the annotation sheet. Sanitize consistently:
Greek α (U+03B1) → `aa`, em/en-dashes → `-`, collapse internal whitespace; RCTD column
names replace non-alphanumerics with `_`. **Don't trust the sheet — fix it in code.**
*Examples from one sheet:* a broad label `"Cd8aa _IEL"` (stray space) silently split into
a second cell type; and `Mid_villus_proximal_enterocytes` carried **two different hex
codes** across its clusters, so first-occurrence dedup picked the wrong color. Define the
palette explicitly in code rather than relying on sheet row order.

**Q: How are clusters mapped to names?**
Via an external CSV: cluster `#` → `Population` (fine, ~41 labels) → `Broad` (~24 labels).
Keep a two-level scheme so you can deconvolve at either granularity without re-annotating.

## Statistical Limitations (read this before claiming significance)

**Q: Why can't my n=2 vs n=2 comparison ever reach significance at the sample level?**
With **2 vs 2** replicates, a two-sided Wilcoxon's **minimum possible p-value is 1/3** —
it *cannot* be significant no matter how large the effect. Meanwhile **spot-level** tests
pool hundreds of thousands of spots and are **pseudoreplicated** (massively inflated
significance). **Neither alone is trustworthy.** *Approach used:* a two-tier rule — call a
result "high-confidence" only if it is **spot-level significant AND consistent in
direction across both replicates** in each group, and corroborate with **pseudobulk
DESeq2** (the correct replication unit). *Example:* Ileum-vs-Duodenum RCTD-proportion and
DEG comparisons.

## Seurat v5 Gotchas

**Q: `AverageExpression` / assay operations behave oddly in Seurat v5.**
v5 returns a **sparse matrix** from `AverageExpression` → wrap in `as.matrix()`. Set
`DefaultAssay()` on the *object you pass* (not a copy). Call `JoinLayers()` before pulling
`counts` from a merged object, or you'll silently get one layer.

## Visualization Conventions

**Q: My composition bars order cell types randomly and colors drift between plots.**
ggplot orders a character fill **alphabetically** by default. Define an **explicit,
ordered palette** (a named color vector) and apply it as **factor levels** to *both* the
fill variable and the dominant-cell-type column so every plot matches. *Example:* the 24
broad types are ordered by compartment (Epithelia → Immune → Stroma, crypt→villus within),
and that same `broad_levels` vector drives legends, stacked bars, and dominant-type maps.

## Bin Size

**Q: Do 8µm and 16µm need different handling?**
Yes. **8µm** applies a `UMI_min = 50` cutoff in RCTD; **16µm** binned data has lower
per-spot UMI, so the cutoff is skipped. They also use different object names and memory
profiles. The config's small-bin default (`BIN_SIZE <- 2`) is **optimistic** — 2µm is far
heavier than the 8/16µm runs this pipeline was exercised on; size resources accordingly.
