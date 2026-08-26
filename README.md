If you use this software in academic work, please cite the associated publication appropriately:

Gaming Addiction Transmission: Interpersonal Neural Pathways Revealed by HYPER-NESS
Chuyang Sun, Mattia Rosso, Ruoyu Niu, Xinyu Ye, Tinghao Tang, Peter Vuust, Leonardo Bonetti, Rixin Tang
bioRxiv 2026.06.18.732549; doi: https://doi.org/10.64898/2026.06.18.732549

# hyperness_fNIRS
Social brain network analysis based on GED (fNIRS), advanced in seperate inter- and intra-brain contribution.

Generalized eigendecomposition for multi-brain (hyperscanning) data: joint covariance across *N* brains, exact decomposition of a component into within-brain and between-brain contributions, and recovery of the corresponding time series.

Three files:

| File | What it does |
|---|---|
| `hyperness_pipeline.m` | the whole flow: data matrix → GED across frequencies → band → decomposition → time series |
| `hyperness_contributions.m` | splits component power into intra- and inter-brain contributions |
| `hyperness_recover.m` | reconstructs the intra- and inter-brain time series of a component |

Works for any number of brains and unequal channel counts: a dyad is `nchan_per_brain = [22 22]`, a triad `[22 22 22]`, and `[22 20 18]` is fine too.

**Dependency:** `filterFGx.m` (Gaussian narrow-band filter, FREQ-NESS external functions) on the MATLAB path. Nothing else.

---

## Method

Stack all channels of all brains into `X` (`nChanTotal × nTime`). For each centre frequency *f*, solve

```
S(f) v = λ Rreg v
```

with `S(f)` the narrow-band covariance, `R` the broadband covariance and `Rreg = (1−γ)R + γ·mean(eig(R))·I`. The eigenvector `v` is a spatial filter defined over all brains at once, so a component can express within-brain activity, between-brain coupling, or a mixture.

**Decomposition.** Writing `v` and `S` in blocks, one per brain,

```
v' S v = Σ_a v_a' S_aa v_a  +  Σ_{a<b} 2 v_a' S_ab v_b
         └──── intra ────┘     └───── inter ─────┘
```

exactly and without remainder. Intra terms are non-negative; inter terms are signed (positive = in phase, negative = anti phase). `hyperness_contributions` verifies the split against `v'*S*v` at every frequency.

**Time series.** Projecting each brain through its own part of the filter and narrow-band filtering gives `s_a(t)`, and

```
s(t)² = Σ_a s_a(t)²  +  Σ_{a<b} 2 s_a(t) s_b(t)
```

is the instantaneous version of the same split: the time average of these series equals the covariance-based ratios exactly. These are the series to feed into Granger causality, cross-frequency coupling, or brain–behaviour analyses.

---

## Usage

```matlab
X = ...;                        % [nChanTotal x nTime], channels x time
nchan_per_brain = [22 22];      % brains stacked in this order
srate = 30.3;

frex   = 0.01:0.002:1;
fwidth = 0.02 * ones(size(frex));
fwidth(frex < 0.02) = 0.01;

% ... GED (section 4 of hyperness_pipeline.m) -> GED struct with
%     .evecs [nChan x nComps x nFreq], .frex, .fwidth, .srate

freq_indices = find(frex >= 0.020 & frex <= 0.026);
comp_idx     = 1;

C  = hyperness_contributions(GED, X, nchan_per_brain, freq_indices, comp_idx);
TS = hyperness_recover(GED, X, nchan_per_brain, freq_indices, comp_idx);
```

`C.intra` and `C.inter` are the band-averaged shares; `C.intra_ratio_f` and `C.inter_ratio_f` give the same quantities frequency by frequency, which is worth checking before trusting a band average. `C.inter_chan{a,b}` is the channel-by-channel inter-brain contribution, whose entries sum to `C.inter(a,b)`.

`TS.intra`, `TS.inter` and `TS.check` hold the time-resolved decomposition and the numerical verification of the identity above.

`hyperness_pipeline.m` runs all of this end to end. Fill in section 1 with your own loader; the script checks that `X`, `nchan_per_brain` and `srate` are defined and stops with a message if they are not. Section 5 marks where the between-condition statistics that define the frequency band would go — they need the whole sample, so they are described in a comment rather than run.

---

## Notes

* **Normalisation.** The decomposition partitions a covariance matrix, so without per-channel z-scoring a brain with larger raw amplitude gets a larger intra-brain share for non-neural reasons. Section 3 of the pipeline does this; set `do_zscore = false` to reproduce the earlier behaviour.
* **Spatial patterns.** Plot `GED.pats` (= `S*v`), not the eigenvector, when interpreting topographies (Haufe et al., 2014, *NeuroImage*, 87, 96–110).
* **Common sources.** Shared physiology (Mayer waves, respiration) and shared stimuli produce inter-brain contributions without any interaction. GED estimates components, it does not assign them a cause; pseudo-pair controls are still needed.
* **Granger causality.** `extras/hyperness_gc.m` implements the directed analysis with a surrogate null that circularly shifts the source series. Do not build the null by shuffling the target series: for narrow-band signals that null collapses to a different scale and declares almost everything significant. See `MIGRATION.md` for the simulation.

---

## References

* Cohen, M. X. (2022). A tutorial on generalized eigendecomposition for denoising, contrast enhancement, and dimension reduction in multichannel electrophysiology. *NeuroImage*, 247, 118809.
* de Cheveigné, A., & Parra, L. C. (2014). Joint decorrelation, a versatile tool for multichannel data analysis. *NeuroImage*, 98, 487–505.
* Rosso, M., Fernández-Rubio, G., Keller, P. E., Brattico, E., Vuust, P., Kringelbach, M. L., & Bonetti, L. (2025). FREQ-NESS reveals the dynamic reconfiguration of frequency-resolved brain networks during auditory stimulation. *Advanced Science*, 2413195.
* Maris, E., & Oostenveld, R. (2007). Nonparametric statistical testing of EEG- and MEG-data. *Journal of Neuroscience Methods*, 164(1), 177–190.
