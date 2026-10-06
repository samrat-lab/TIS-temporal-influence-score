# Input data

Both input files are **included in this repository**, so the pipeline runs
out of the box:

```
data/
  data_info.mat         1.5 MB   variable: data_info
  Final_TF_Target.mat    23 MB   variable: Final_TF_Target
```

From the repository root:

```matlab
out = run_all();
```

No paths need to be edited. `config/tis_config.m` resolves both files
relative to this folder.

---

## 1. `data_info.mat` - the proteomics measurements

`data_info` is a `{12 x 2}` cell array, **one row per post-infection
timepoint**, in the order

```
6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50 hours
```

```
data_info                       {12 x 2} cell
 |
 +-- data_info{t,1}             timepoint label
 |
 +-- data_info{t,2}             {4 x >=6} cell  <- one ROW per strain
                                                   one COLUMN per protein list
```

### Strain rows of `data_info{t,2}`

| row | strain | note |
|---|---|---|
| 1 | H37Ra  | attenuated laboratory strain |
| 2 | H37Rv  | virulent laboratory strain |
| 3 | BND433 | clinical isolate |
| 4 | JAL2287| clinical isolate |

This row order must match `cfg.strain_labels` in `config/tis_config.m`.

### Columns of `data_info{t,2}`

Each cell is a string column vector of protein identifiers. Columns beyond
the sixth, if present, are ignored by the pipeline.

| col | contents | used by the pipeline |
|---|---|---|
| 1 | master protein list | no |
| 2 | **newly synthesized proteins detected** at this timepoint in this strain | **yes** - defines `V_t`, the node set |
| 3 | up-regulated proteins (infected vs uninfected control) | no |
| 4 | down-regulated proteins (infected vs uninfected control) | no |
| 5 | **DENSP** = union of columns 3 and 4 | optional, see below |
| 6 | **transcription factors detected** at this timepoint | **yes** - defines `T_t`, the TF subset |

Fold changes underlying columns 3-5 are computed for each strain against the
uninfected control at the same timepoint. See the Methods section of the
paper for the experimental protocol and the data-processing rules.

### How the pipeline uses these columns

Only two columns enter the score:

```
V_t = data_info{t,2}{k,2}          all proteins measured at timepoint t
T_t = data_info{t,2}{k,6}          transcription factors among them
P_t = V_t \ T_t                    newly synthesized NON-TF proteins
```

`P_t` is the denominator of the local influence term `L`, and `T_t` supplies
the nodes over which the propagated term `N` is summed.

To score on DENSP instead of the full newly synthesized proteome, change one
line in `config/tis_config.m`:

```matlab
cfg.col_protein_pool = cfg.col_densp;      % use column 5 instead of column 2
```

or override it at the call site: `run_all('col_protein_pool', 5)`.

### Inspecting the file

```matlab
load('data/data_info.mat')

size(data_info)                    % 12 x 2  (timepoints)
size(data_info{1,2})               % 4 x ... (strains x columns)

% newly synthesized proteins, H37Rv (row 2), 6 h (row 1)
p = data_info{1,2}{2,2}(:,1);
numel(p)
p(1:5)

% transcription factors in the same sample
tf = data_info{1,2}{2,6}(:,1);
numel(tf)
tf(1:5)

% protein counts per timepoint for every strain
cellfun(@(blk) numel(blk{2,2}(:,1)), data_info(:,2))'
```

---

## 2. `Final_TF_Target.mat` - the regulatory backbone

`Final_TF_Target` is an `[n x 2]` string array of curated transcription
factor to target interactions taken from the **TFLink** database.

| column | contents |
|---|---|
| 1 | source: transcription factor |
| 2 | target: regulated gene / protein |

```matlab
load('data/Final_TF_Target.mat')
size(Final_TF_Target)
Final_TF_Target(1:5,:)
```

The edge list is **time-invariant**: it states which regulatory interactions
are possible, not when they occur. All temporal structure in the analysis
comes from the measurements - an edge is realised at transition `i -> i+1`
only if the source TF was detected at timepoint `i` and the target was
detected at timepoint `i+1`.

The pipeline first restricts this list, per strain, to edges whose source is
a TF detected at some point in that strain (`build_tf_target_info.m`); the
rest can never contribute and are dropped.

Please cite TFLink alongside this work if you reuse the edge list.

---

## 3. Identifier conventions

All identifiers are compared case-insensitively. `src/normalize_ids.m`
uppercases and trims every name before any set operation, so `rv0001`,
`Rv0001` and `RV0001 ` are one entity. The two files must use the same
identifier namespace, otherwise the intersection between the measured
proteome and the TFLink edge list will be empty and every score will be zero.

A quick sanity check that the namespaces line up:

```matlab
load('data/data_info.mat'); load('data/Final_TF_Target.mat');
prot = upper(data_info{1,2}{2,2}(:,1));          % H37Rv, 6 h
sum(ismember(upper(Final_TF_Target(:,2)), prot)) % should be >> 0
```

---

## 4. Replacing the data

To run the same analysis on a different experiment, build a `data_info` cell
with the layout above and update `config/tis_config.m`:

* `strain_labels` - one label per row of `data_info{t,2}`, in row order
* `tp_hours` - one entry per row of `data_info`, in row order
* `col_protein_pool` / `col_tf` - if your columns sit elsewhere

Nothing else in the code assumes four strains or twelve timepoints;
`cfg.n_strain` and `cfg.n_tp` are derived from those two fields.
