# Agent design: how omicos was used in this analysis

This document describes how omicos was used as the analysis assistant
for the gsMap AD whole-body spatial project that this repository organizes. It is
deliberately mechanical — what tools, what patterns, what worked, what didn't —
so other researchers can adopt the same setup.

---

## 1. Agent setup

- **omicos CLI** as the primary interface. Model selection was
  task-dependent: Sonnet for routine edits (notebook patches, file moves,
  inventory scripts) and Opus for multi-step reasoning (figure methodology
  audits, story drafting, causal-evidence synthesis).
- **Stateless code agents with persistent file-based memory** at
  `.omicos/memory/`. Each invocation starts clean; continuity comes from
  on-disk notes (project summary, feedback log, references), not from a
  long-running session.
- **Sub-agents** were used to parallelize work that didn't need shared state:
  - `Explore` sub-agents to scan directory trees, read many files, and report
    back inventories (e.g., "list every gsMap output under `results/`").
  - `Plan` sub-agents to draft multi-step execution plans before any file
    was touched.
  - `general-purpose` sub-agents to fan out independent research tasks
    (literature lookups, schema checks) and return short structured answers.

---

## 2. Interaction patterns

The collaboration converged on a small set of repeated patterns:

- **Iterative refinement.** Small commits, dry-run a cell in `/tmp/cell.py`
  before writing it into the notebook, then render and check the resulting
  PNG. No "big bang" rewrites.
- **Read the PNG, do not just look at the code.** omicos was explicitly
  asked to `Read` rendered figure PNGs and describe what it saw, rather
  than predicting the figure from the source. This caught at least one
  confabulated description (see Section 5).
- **`AskUserQuestion` at branch points.** Whenever 2-3 reasonable paths
  existed (e.g., "merge h5ad lightweight, medium, or full expression
  matrix?"), omicos posed an explicit multiple-choice question instead
  of guessing.
- **Cell IDs for surgical notebook updates.** Edits were targeted by
  stable `cell_id`, looked up via
  `next(i for i,c in enumerate(nb['cells']) if c.get('id')==target)`.
  No reindexing, no rebuilding the notebook.
- **Sanity-check by re-loading JSON after every notebook write.** After
  patching a notebook, immediately `json.load` it again and verify the
  expected source text is present. Several "but I just saved it" bugs
  came from skipping this step.

---

## 3. Memory and skills

### Memory layout

`.omicos/memory/` is organized into four categories:

- **User memory** — durable facts about the operator (email, working
  directory, model preferences).
- **Feedback memory** — short notes capturing corrections from the user
  ("don't rebuild the whole notebook", "use union not intersection for
  organs") so the same mistake isn't repeated next session.
- **Project memory** — one entry per analysis. The relevant entry here
  is **gsMap AD whole-body spatial**: 126 gsMap results across
  15 organs × 10 ages, two-pathway discovery (myeloid + gut epithelial
  APP processing), multi-trait comparison running.
- **Reference memory** — external pointers (paper IDs, dataset GEO
  accessions, tool docs) that recur across projects.

### Skills

Three skills were active for this work:

- `omicos-guide` — house style for omicos interactions in this
  lab (tool-call conventions, file-path discipline).
- `paper-writer` — IMRAD scaffolding, literature matrix, reference
  management. Used to bootstrap `paper/AD_colon_story.md`.
- `omicos-api` — for any auxiliary scripts that themselves call the
  omicos team SDK (e.g., batch description generation).

---

## 4. Workflows

Multi-step work was driven by the **Workflow tool**, which fans out
parallel sub-agents and then joins their outputs. Workflows are
declarative: each step lists its inputs, the agent to invoke, and what
the next step expects.

The exact workflow that produced *this repository* is
`organize-mouse-spatial-repo`:

1. **Inventory** — Explore sub-agent walks the source analysis tree and
   emits a JSON manifest of figures, notebooks, scripts, results, and
   data files.
2. **Classify** — general-purpose sub-agent tags each artifact as
   `keep / archive / regenerate` based on a rubric.
3. **Plan layout** — Plan sub-agent proposes the target directory
   structure (`figures/`, `notebooks/`, `results/`, `scripts/`,
   `data_manifest/`, `docs/`).
4. **Execute moves** — code agent performs the file moves with
   `git mv` where applicable and rewrites notebook paths.
5. **Generate docs** — code agent writes per-directory READMEs and
   this AGENT_DESIGN.md, sourcing context from
   `paper/conversation_walkthrough.md`.
6. **Verify** — code agent reopens each notebook and confirms relative
   paths still resolve.

Steps 1-2 run in parallel; 3-6 run sequentially.

---

## 5. What worked

- **Reality-check rendered PNGs prevented confabulation.** The original
  Fig 4A narrative was "APP-pathway genes colon-strong, brain-absent".
  Asking omicos to actually `Read` the PNG revealed that Whole_Blood
  often beat Colon (APP 13.1, ADAM10 20.5) and that SORL1/TREM2/SPI1
  were uniformly zero in GTEx bulk eQTL. The corrected description
  ("peripheral-strong, brain-absent") is now what the paper says.
- **Strict cell-ID-based updates prevented notebook rebuilds.**
  Targeting `cell_id` and editing only `source` kept notebook diffs
  small and reviewable, and avoided losing outputs.
- **Pre-flight inventory then plan-then-execute reduced waste.** A
  short Explore pass before any file move surfaced naming
  inconsistencies and let the Plan step propose a layout that needed
  no rework.
- **`AskUserQuestion` at the h5ad-merge branch** ("lightweight metadata
  only / medium / full expression + metadata") let the user pick the
  5.57 GB option deliberately rather than discovering it after the run.
- **Methodology audits as two-dimensional checks** — "data correctness"
  and "journal visual conventions" scored separately — caught the
  Fig 4B (PP.H4-only bar instead of stacked H0-H4), 4C (missing coloc
  vs SharePro pairing), and 4D (SMR with no HEIDI) issues.

---

## 6. What didn't work (be honest)

- **Initial size legend in pies didn't match physical scale.** The
  first marsilea pie-on-heatmap version drew legend circles in data
  units while the main pies were in inches. Legend looked plausible
  but was wrong by ~2x. Required a refactor to enforce
  `size_ax_width_inches = 3 * CELL_inches`, `xlim=(0, 3)`,
  `aspect='equal'`. Lesson: when a legend encodes size, share the
  physical (inch) unit with the main artist, not just the data unit.
- **Organ name mismatches bit us several times.** `Brown_Fat` vs `BAT`,
  `Small_Intestine` vs `Ileum`, intersection-vs-union for the
  Young/Old organ set. Each surfaced as "the two figures don't show
  the same organs" and required a mapping dict and explicit
  symmetric-difference diagnostics
  (`set(A) ^ set(B)`) to fix.
- **Some patches failed silently.** Unicode quote characters
  (curly `'` vs straight `'`) in `old_string` parameters meant Edit
  calls reported success but matched nothing. Caught only by the
  "re-load JSON and grep" sanity step. Easy to miss without that
  habit.
- **Cross-dataset assumptions slipped through.** omicos once asserted
  Fig 2A used raw 6w ST when it actually used CMAP multi-age. Required
  the user to push back, and omicos to verify by inspecting the
  `sample` column. Turned into a useful methodological point (raw ST
  spots dilute epithelial signal, so CMAP projection is necessary),
  but the initial confidence was unwarranted.

---

## 7. Reusability

Other researchers can adopt the same setup with three concrete steps:

1. **Build the memory directory first.** Before any analysis, create
   `.omicos/memory/` with at least a one-paragraph project summary
   and a feedback log. Update both as you go. The agent is stateless;
   the memory is the project.
2. **Use `AskUserQuestion` at branch points.** Whenever you find
   yourself about to guess between 2-3 reasonable defaults
   (file size tradeoff, intersection vs union, which organ set,
   which colormap), instruct the agent to ask. The cost of a
   question is far less than the cost of redoing a 5 GB merge.
3. **Render-check at every visualization step.** Have the agent
   read the rendered PNG and describe what it sees before you write
   any caption or interpretation. This single habit prevents the
   most common failure mode (plausible-sounding figure descriptions
   that don't match the actual image).

The same three habits — memory, branch-point questions, render-checks —
generalize beyond spatial transcriptomics to any analysis with
iterative figure refinement.
