# Temporal Influence Score (TIS)

MATLAB implementation of the temporal influence score used to rank
transcription factors in time-resolved infection proteomics. The pipeline
builds a directed, layered TF -> newly-synthesized-protein network across
post-infection timepoints, computes the two components of the score
recursively, and tests them against a permutation null.

For a transcription factor `t` detected at timepoint `i`,

```
S(i)\_t  =  L(i->i+1)\_t  +  N(i->i+1)\_t

           |N\_target(t)|          sum\_{u in N\_TF(t)} S(i+1)\_u
        =  --------------   +    -----------------------------
             |P\_{i+1}|           sum\_{v in T\_{i+1}} S(i+1)\_v

           1st term: local         2nd term: network-propagated
```

where `P\_{i+1}` is the set of newly synthesized non-TF proteins at the next
timepoint, `N\_target(t)` those of them that `t` regulates, `T\_{i+1}` the
transcription factors at the next timepoint and `N\_TF(t)` those of them that
`t` regulates. The recursion runs backwards in time from the terminal
transition, where no downstream TF layer exists and `S = L`.

Full equation-to-code mapping: [`docs/METHODS.md`](docs/METHODS.md).

## Requirements

* MATLAB R2019b or newer (string arrays, `digraph`, `writetable`)
* No toolboxes beyond base MATLAB

## Quick start

**The input data ships with this repository** (`data/data\_info.mat` and
`data/Final\_TF\_Target.mat`), so nothing needs to be downloaded or
reconfigured. From the repository root:

```matlab
out = run\_all();                              % full analysis, 1000 permutations

out = run\_all('n\_perm', 100);                 % faster smoke run
out = run\_all('run\_permutation', false);      % scores only
```

To run the same pipeline on data held elsewhere:

```matlab
out = run\_all('data\_dir', '/path/to/folder');                 % both files together
out = run\_all('data\_info\_file', '/path/to/data\_info.mat', ... % files apart
              'tflink\_file',    '/path/to/Final\_TF\_Target.mat');
```

Check the implementation against the hand-computed toy network first:

```matlab
addpath('tests'); test\_tis\_toy
```

## Repository layout

```
run\_all.m                        end-to-end driver
config/tis\_config.m              paths, strain order, timepoints, column map
src/
  load\_inputs.m                  load and validate the two .mat inputs
  normalize\_ids.m                canonical identifier form
  tag\_nodes.m / untag\_nodes.m    "NAME\_<timepoint>" node ids, for net.edges
  build\_tf\_target\_info.m         strain TF universe + TFLink sub-network
  build\_node\_sets.m              V\_t, T\_t, P\_t per timepoint
  compute\_tis\_scores.m           recursive L, N, S          <- core of the method
  build\_temporal\_network.m       layered directed edge set (Section 2.3)
  assemble\_score\_matrices.m      TF x timepoint score matrices
  run\_permutation\_test.m         empirical p-values
  export\_results.m               CSV tables
  export\_legacy\_mats.m           original cell-array .mat layout
tests/test\_tis\_toy.m             unit test with hand-computed expectations
docs/METHODS.md                  equations -> code
data/
  data\_info.mat                  proteomics measurements (included)
  Final\_TF\_Target.mat            TFLink TF -> target edge list (included)
  README.md                      full description of both files
```

## Configuration

Everything that depends on the experiment lives in `config/tis\_config.m`:

|field|default|meaning|
|-|-|-|
|`strain\_labels`|`{'H37Ra','H37Rv','BND433','JAL2287'}`|row order of `data\_info{t,2}`|
|`tp\_hours`|`\[6 10 ... 50]`|post-infection timepoints (h)|
|`col\_protein\_pool`|`2`|`data\_info` column defining `V\_t`; set to `5` to score on DENSP|
|`col\_tf`|`6`|`data\_info` column listing detected TFs|
|`n\_perm`|`1000`|permutations|
|`alpha`|`0.05`|significance threshold|
|`rng\_seed`|`42`|random seed|

Any field can be overridden at the call site, e.g.
`run\_all('col\_protein\_pool', 5)`.

## Inputs

Two `.mat` files, both included in `./data` and documented in detail in
[`data/README.md`](data/README.md):

**`data\_info.mat`** - the proteomics measurements. A `{12 x 2}` cell, one row
per post-infection timepoint (6, 10, ... 50 h). `data\_info{t,2}` is a
`{4 x >=6}` cell with **one row per strain** (H37Ra, H37Rv, BND433, JAL2287)
and one column per protein list:

|col|contents|used|
|-|-|-|
|1|master protein list||
|2|newly synthesized proteins detected at this timepoint in this strain|defines `V\_t`|
|3|up-regulated proteins (vs uninfected control)||
|4|down-regulated proteins (vs uninfected control)||
|5|DENSP = union of columns 3 and 4|optional node set|
|6|transcription factors detected at this timepoint|defines `T\_t`|

Fold changes are computed per strain against the uninfected control at the
same timepoint. The score uses columns 2 and 6; the non-TF protein layer is
`P\_t = V\_t \\ T\_t`.

**`Final\_TF\_Target.mat`** - the regulatory backbone. An `\[n x 2]` string
array of TFLink interactions, column 1 the source transcription factor and
column 2 its target. The list is time-invariant: it says which interactions
are possible, not when they occur. An edge is realised at transition
`i -> i+1` only if the source was detected at `i` and the target at `i+1`.

Identifiers are compared case-insensitively; both files must use the same
namespace. `data/README.md` gives inspection snippets and a namespace sanity
check.

## Outputs

```
results/
  tis\_pipeline.mat                       complete workspace
  TIS\_summary.csv                        TFs and significant TFs per strain
  scores/
    TIS\_scores\_<strain>.csv              TF x (L, N, S) for all 12 timepoints
  permutation/
    TIS\_pvalues\_<strain>.csv             observed scores, p-values, null means
    TIS\_significant\_TFs\_<strain>.csv     TFs with p < alpha at >= 1 transition
  legacy/
    score\_mat\_1st\_term.mat               original cell-array layout
    score\_mat\_2nd\_term.mat
    p\_mat\_1st\_term.mat
```

`TIS\_scores\_<strain>.csv` has one row per transcription factor and three
blocks of twelve columns, `L\_6h ... L\_50h`, `N\_6h ... N\_50h` and
`S\_6h ... S\_50h`. The 50 h columns are zero by construction: no transition
starts at the last timepoint.

The directed temporal network itself is returned in
`out.net{k}.edges` as a tagged `\[E x 2]` edge list (`NAME\_<timepoint>`), for
inspection or export to a network viewer.

## 

