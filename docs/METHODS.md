# Methods: from the equations to the code

This note maps every equation of the Methods section onto the function that
implements it. Symbols follow the paper.

| symbol | meaning | where it lives in the code |
|---|---|---|
| `i` | transition index, `1 ... n_tp-1` | loop variable in `compute_tis_scores.m` |
| `V_i` | all proteins measured at timepoint `i` | `nodes{k}.prot{i}` |
| `T_i` | transcription factors detected at timepoint `i` | `nodes{k}.tf{i}` |
| `P_i` | newly synthesized non-TF proteins, `V_i \ T_i` | `nodes{k}.nonTF{i}` |
| `E_i` | temporal edges `i -> i+1` | `net{k}.layer_edges{i}` |
| `L^(i->i+1)_t` | local influence (1st term) | `tis{k}.L{i}`, `scores{k}.L` |
| `N^(i->i+1)_t` | network-propagated influence (2nd term) | `tis{k}.N{i}`, `scores{k}.N` |
| `S^(i)_t` | temporal influence score | `tis{k}.S{i}`, `scores{k}.S` |
| `p^(i)_t` | empirical p-value | `perm{k}.p_S` |

## 1. Node set (Eq. 1)

`V_i = T_i u P_i`. Implemented in `build_node_sets.m`. The protein pool is
column 2 of `data_info` (newly synthesized proteome); the TF subset is column
6. `P_i` is obtained by set difference, so a protein listed in both columns is
counted once, as a TF.

## 2. Temporal edges (Eqs. 2-4)

`build_temporal_network.m`. Two edge types are formed from the time-invariant
TFLink list:

* TF -> non-TF: `E_{t->p} = {(t,p) : t in T_i, p in P_{i+1}, (t->p) in TFLink}`
* TF -> TF:     `E_{t->u} = {(t,u) : t in T_i, u in T_{i+1}, (t->u) in TFLink}`

An edge is realised only if BOTH endpoints were measured, the source at layer
`i` and the target at layer `i+1`. The active source set at layer `i` is the
union of TFs reached from layer `i-1` and TFs newly detected at layer `i`.

Nodes are `(protein, timepoint)` pairs written `NAME_i` (`tag_nodes.m`), so
the same protein at two timepoints is two nodes and the graph is acyclic by
construction.

## 3. Neighbourhoods (Eqs. 5-7)

For a TF `t` at layer `i`, the forward neighbourhood `N_t` is its set of
TFLink targets measured at `i+1`. It is partitioned into

* `N_target(t) = N_t n P_{i+1}` - non-TF targets,
* `N_TF(t)     = N_t n T_{i+1}` - downstream TF targets.

Both are computed inside the TF loop of `compute_tis_scores.m` with
`intersect` / `setdiff` on normalized identifiers, so each distinct target is
counted once.

## 4. Local influence, 1st term (Eq. 10)

```
L^(i->i+1)_t = |N_target(t)| / |P_{i+1}|
```

The denominator is the complete non-TF layer, identical for every TF at that
transition, so `L` is directly comparable across TFs within a timepoint. If
`|P_{i+1}| = 0` the term is defined as 0.

## 5. Network-propagated influence, 2nd term (Eq. 11)

```
N^(i->i+1)_t = sum_{u in N_TF(t)} S^(i+1)_u / sum_{v in T_{i+1}} S^(i+1)_v
```

The denominator runs over all TFs at layer `i+1` that have a defined score.
Because the term depends on scores one layer ahead, the recursion is
evaluated backwards, `i = n_tp-1` down to `1`.

## 6. Base case (Eq. 9)

At the terminal transition `i = n_tp-1` there is no layer `i+2` from which
downstream TF scores could be drawn, so `N = 0` and `S = L`.

## 7. Total score (Eqs. 8, 12)

`S^(i)_t = L^(i->i+1)_t + N^(i->i+1)_t`. Reported as `scores{k}.S`, a
`nTF x n_tp` matrix per strain. The final column is zero throughout: no
transition starts at the last timepoint.

`assemble_score_matrices.m` writes these values directly from the per-
timepoint vectors returned by the recursion. No filtering against the network
node list is applied, and none is needed: a TF that is absent at a timepoint,
or whose TFLink targets were all unmeasured at the next one, has
`L = N = 0` by definition, which is exactly the zero the matrix is
initialised with.

## 8. Permutation test (Eq. 13)

`run_permutation_test.m`. For each strain and each transition, the observed
node set at layer `i+1` is replaced by an equally sized sample drawn without
replacement from the union of proteins detected across all timepoints in that
strain. The TF -> target topology is untouched. `L` and `N` are recomputed and

```
p^(i)_t = (1/n) sum_k 1[ S^(i)_{t,null,k} >= S^(i)_t ],   n = 1000
```

TFs with `p < 0.05` at one or more transitions are retained.

Two properties of this null are deliberate and should be stated when the
result is interpreted:

1. A randomly drawn protein has no timepoint of its own, so its TF status is
   taken from the strain-wide TF universe rather than from the per-timepoint
   TF list.
2. The propagated term re-uses the **observed** downstream scores as weights;
   only the membership of the next layer is randomised. The test therefore
   asks whether a TF reaches an unusually influential part of the next layer,
   not whether the entire recursive cascade could arise by chance.

p-values are computed as `count / n_perm`; a value of exactly 0 is reported as
`1/n_perm` in the exported tables, the resolution limit of the test.

### Implementation note

The inner loop is expressed as two sparse matrix-vector products. With
`M_non` and `M_tf` the `nTF x n_pool` indicator matrices of each TF's non-TF
and TF targets, and `x` the 0/1 membership vector of one permuted layer,

```
L_null = (M_non * x) / sum(x & ~is_TF)
N_null = (M_tf  * (x .* w)) / sum(w)
```

where `w` holds the observed scores of layer `i+1` indexed over the pool.
This is algebraically identical to looping over TFs and is what makes 1000
permutations tractable.

## Relationship to Supplementary Note 1

The worked example in the supplement specifies a separate edge list for each
transition, which makes the arithmetic easy to follow by hand. The pipeline
instead holds one time-invariant TFLink edge list and derives the per-
transition edges by intersecting it with what was measured at each timepoint.
The unit test in `tests/test_tis_toy.m` therefore uses its own small network,
built so that the measured sets alone determine the realised edges, and checks
the same three cases the supplement illustrates: the terminal transition, a
purely local transition, and a transition combining local and propagated
influence.
