%% ====================================================================
%  HyperNESS pipeline: multi-brain data -> GED -> intra/inter -> time series
%
%  Works for any number of brains: a dyad is nchan_per_brain = [22 22],
%  a triad is [22 22 22]. Channel counts may differ between brains.
%
%  Dependency: filterFGx.m (FREQ-NESS external functions) on the MATLAB path.
% ====================================================================

clear; clc;

addpath(fileparts(mfilename('fullpath')));
% addpath('D:\Matlab\toolbox\FREQNESS_Toolbox\FREQNESS_ExternalFunctions\');

%% ====================================================================
%  1. DATA
%
%  Fill in this section. Three variables are needed:
%
%    X               [nChanTotal x nTime], channels x time (note the
%                    orientation), brains stacked in the order given by
%                    nchan_per_brain, already segmented to the condition of
%                    interest - the covariance matrices are estimated on
%                    whatever you pass in
%    nchan_per_brain [22 22] for a dyad, [22 22 22] for a triad; channel
%                    counts may differ between brains
%    srate           sampling rate in Hz
%
%  A typical loader:
%
%    raw   = load(fullfile(dataPath, 'dyad07.mat'));
%    seg   = raw.data(raw.events.task1Start:raw.events.task1End, :)';   % 44 x T
%    X     = seg;
%    nchan_per_brain = [22 22];
%    srate = 30.3;
%
%  Interpolate or reject NaNs before this point.
% ====================================================================

% <- your loading code here

if ~exist('X', 'var') || ~exist('nchan_per_brain', 'var') || ~exist('srate', 'var')
    error('Section 1: define X, nchan_per_brain and srate before running the script.');
end
if sum(nchan_per_brain) ~= size(X, 1)
    error('sum(nchan_per_brain) = %d does not match size(X,1) = %d.', ...
        sum(nchan_per_brain), size(X, 1));
end
if any(isnan(X(:)))
    error('X contains NaNs; interpolate or reject them first.');
end

fprintf('%d brains, %d channels, %d samples (%.1f min)\n', ...
    numel(nchan_per_brain), size(X, 1), size(X, 2), size(X, 2) / srate / 60);

%% ====================================================================
%  2. PARAMETERS
% ====================================================================

frex   = 0.01:0.002:1;                  % frequency grid
fwidth = 0.02 * ones(size(frex));       % Gaussian FWHM per frequency
fwidth(frex < 0.02) = 0.01;             % narrower below 0.02 Hz

ncomps         = 10;
regularisation = 0.01;

%% ====================================================================
%  3. NORMALISATION
%  Per-channel z-score. The decomposition partitions a covariance matrix, so
%  without this a brain with larger raw amplitude gets a larger intra share
%  for purely non-neural reasons. Set to false to reproduce the old pipeline.
% ====================================================================

do_zscore = true;

X = X - mean(X, 2);
if do_zscore
    sd = std(X, 0, 2);
    sd(sd == 0 | ~isfinite(sd)) = 1;
    X = X ./ sd;
end

%% ====================================================================
%  4. GED, one decomposition per frequency
%  S(f) v = lambda * Rreg v, with S the narrow-band covariance and R the
%  broadband covariance. Replace this block by
%
%     GED = FREQNESS_NetworkEstimation(X, frex, srate, 'regularisation', ...
%              regularisation, 'ncomps', ncomps, 'filter', 'logarithmic', ...
%              'fwidth', fwidth);
%     GED.fwidth = fwidth;
%
%  if you prefer the FREQ-NESS implementation. The two functions below only
%  need GED.evecs, GED.frex, GED.fwidth and GED.srate.
% ====================================================================

nchan = size(X, 1);

covR = cov(X');
covR = (covR + covR') / 2;
covR = (1 - regularisation) * covR + regularisation * mean(eig(covR)) * eye(nchan);

% whitening transform: solving the symmetric problem in the whitened space is
% more stable than eig(covS, covR) and keeps eigenvalues real
[Vr, Dr] = eig(covR);
[dr, o]  = sort(diag(Dr), 'descend');
Vr = Vr(:, o);
keep = dr > max(dr) * 1e-10;
Wht  = Vr(:, keep) * diag(1 ./ sqrt(dr(keep)));

ncomps = min(ncomps, size(Wht, 2));
evecs  = zeros(nchan, ncomps, numel(frex));
evals  = zeros(ncomps, numel(frex));
pats   = zeros(nchan, ncomps, numel(frex));

fprintf('GED over %d frequencies ', numel(frex));
for fi = 1:numel(frex)

    narrowData = filterFGx(X, srate, frex(fi), fwidth(fi), 0);
    covS = cov(narrowData');
    covS = (covS + covS') / 2;

    Sw = Wht' * covS * Wht;
    Sw = (Sw + Sw') / 2;
    [Vs, Ds] = eig(Sw);
    [lam, o] = sort(diag(Ds), 'descend');

    V = Wht * Vs(:, o(1:ncomps));
    V = V ./ sqrt(sum(V.^2, 1));            % unit-norm spatial filters

    evecs(:, :, fi) = V;
    evals(:, fi)    = lam(1:ncomps);
    pats(:, :, fi)  = covS * V;             % activation patterns, for topographies

    if mod(fi, 50) == 0, fprintf('.'); end
end
fprintf(' done\n');

GED = struct('evecs', evecs, 'evals', evals, 'pats', pats, ...
             'frex', frex, 'fwidth', fwidth, 'srate', srate);

%% ====================================================================
%  5. FREQUENCY BAND AND COMPONENT
%
%  NOT RUN HERE. In a real study the band comes from a contrast between
%  conditions: loop sections 1-4 over every dyad and condition, collect
%  GED.evals into a [dyads x frequencies] matrix per condition, contrast the
%  conditions frequency by frequency and correct over frequencies with a
%  cluster-based permutation test (Maris & Oostenveld, 2007, J Neurosci
%  Methods, 164, 177-190). The surviving cluster defines freq_range.
%
%  Here a band is simply declared so that sections 6-7 can run.
% ====================================================================

freq_range   = [0.020 0.026];
comp_idx     = 1;
freq_indices = find(frex >= freq_range(1) & frex <= freq_range(2));

fprintf('band %.3f-%.3f Hz (%d frequency points), component %d\n', ...
    freq_range(1), freq_range(2), numel(freq_indices), comp_idx);

%% ====================================================================
%  6. INTRA / INTER DECOMPOSITION
% ====================================================================

C = hyperness_contributions(GED, X, nchan_per_brain, freq_indices, comp_idx);

for a = 1:numel(nchan_per_brain)
    fprintf('intra brain %d : %6.1f %%\n', a, 100 * C.intra(a));
end
for a = 1:numel(nchan_per_brain)
    for b = (a+1):numel(nchan_per_brain)
        fprintf('inter %d-%d     : %6.1f %%\n', a, b, 100 * C.inter(a, b));
    end
end
fprintf('total         : %6.1f %%\n', 100 * (C.intra_total + C.inter_total));

% C.intra_ratio_f / C.inter_ratio_f hold the same quantities frequency by
% frequency: check them before trusting the band average, since a sign flip
% inside the band makes the average meaningless.
% C.inter_chan{1,2} is the channel-by-channel inter-brain contribution.

%% ====================================================================
%  7. TIME-SERIES RECOVERY
% ====================================================================

TS = hyperness_recover(GED, X, nchan_per_brain, freq_indices, comp_idx);

fprintf('instantaneous identity, max deviation: %.3g\n', max(TS.check));

% the time average of the recovered series reproduces section 6
k = 1;
Ptot = mean(TS.component(1, :, k).^2);
fprintf('at %.3f Hz: intra brain 1 = %.4f (time series) vs %.4f (covariance)\n', ...
    TS.freqs(k), mean(TS.intra(1, :, k)) / Ptot, C.intra_ratio_f(1, k));
fprintf('at %.3f Hz: inter 1-2     = %.4f (time series) vs %.4f (covariance)\n', ...
    TS.freqs(k), mean(TS.inter(1, :, k)) / Ptot, C.inter_ratio_f(1, 2, k));

%% ====================================================================
%  8. SAVE
% ====================================================================

% save('hyperness_results.mat', 'GED', 'C', 'TS', 'freq_range', ...
%      'comp_idx', 'nchan_per_brain', 'frex', 'fwidth', 'srate', '-v7.3');
