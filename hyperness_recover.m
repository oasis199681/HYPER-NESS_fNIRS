function TS = hyperness_recover(GED, X, nchan_per_brain, freq_indices, comp_idx)
% HYPERNESS_RECOVER  Reconstruct the intra- and inter-brain time series of a
% multi-brain GED component.
%
%   TS = hyperness_recover(GED, X, nchan_per_brain, freq_indices, comp_idx)
%
% INPUT   same as HYPERNESS_CONTRIBUTIONS. freq_indices may be a vector, in
%         which case every output gets a third dimension, one slice per
%         frequency.
%
% METHOD
%   1. project the broadband data of each brain through its own part of the
%      spatial filter:            s_a(t) = v_a' * X_a(t)
%   2. narrow-band filter each projection (projection and filtering are both
%      linear, so the order does not matter)
%   3. split the instantaneous component power
%
%        s(t) = SUM_a s_a(t)
%        s(t)^2 = SUM_a s_a(t)^2 + SUM_{a<b} 2 s_a(t) s_b(t)
%                 \___ intra ___/  \_______ inter _______/
%
%   This is the time-resolved version of HYPERNESS_CONTRIBUTIONS: the time
%   average of these series equals the covariance-based ratios exactly.
%   (The original dyadic code used s1.*s2 without the factor 2, which only
%   rescales the inter series.)
%
% OUTPUT   nB brains, nP = nB*(nB-1)/2 pairs, T samples, nF frequencies
%   TS.s           [nB x T x nF]  narrow-band per-brain component signals
%   TS.component   [1  x T x nF]  full component signal, sum over brains
%   TS.intra       [nB x T x nF]  per-brain instantaneous power
%   TS.inter       [nP x T x nF]  per-pair instantaneous cross term
%   TS.inter_total [1  x T x nF]  sum over pairs
%   TS.pairs       [nP x 2]       brain indices of each pair
%   TS.check       [1 x nF]       max |component^2 - (intra + inter)|, ~1e-13
%   TS.freqs, TS.comp, TS.srate
%
% DEPENDENCY: filterFGx.m (FREQ-NESS external functions)

nBrain = numel(nchan_per_brain);
nFreq  = numel(freq_indices);
T      = size(X, 2);

idx = cell(1, nBrain);
o = 0;
for a = 1:nBrain
    idx{a} = (o + 1):(o + nchan_per_brain(a));
    o = o + nchan_per_brain(a);
end
if o ~= size(X, 1)
    error('sum(nchan_per_brain) = %d does not match size(X,1) = %d.', o, size(X, 1));
end

pairs = zeros(nBrain*(nBrain-1)/2, 2);
q = 0;
for a = 1:nBrain
    for b = (a+1):nBrain
        q = q + 1;
        pairs(q, :) = [a b];
    end
end
nPair = size(pairs, 1);

s          = zeros(nBrain, T, nFreq);
intra      = zeros(nBrain, T, nFreq);
inter      = zeros(nPair,  T, nFreq);
component  = zeros(1, T, nFreq);
check      = zeros(1, nFreq);

for k = 1:nFreq

    freq_idx = freq_indices(k);
    v = GED.evecs(:, comp_idx, freq_idx);

    % ---- broadband projection, one per brain ----
    s_broad = zeros(nBrain, T);
    for a = 1:nBrain
        s_broad(a, :) = real(v(idx{a})' * X(idx{a}, :));
    end

    % ---- narrow-band filtering, same width as the GED ----
    s_filt = filterFGx(s_broad, GED.srate, GED.frex(freq_idx), GED.fwidth(freq_idx), 0);
    s_filt = s_filt - mean(s_filt, 2);

    s(:, :, k) = s_filt;
    component(1, :, k) = sum(s_filt, 1);

    % ---- instantaneous power decomposition ----
    for a = 1:nBrain
        intra(a, :, k) = s_filt(a, :).^2;
    end
    for q = 1:nPair
        a = pairs(q, 1);
        b = pairs(q, 2);
        inter(q, :, k) = 2 * s_filt(a, :) .* s_filt(b, :);
    end

    check(k) = max(abs(component(1, :, k).^2 - ...
        (sum(intra(:, :, k), 1) + sum(inter(:, :, k), 1))));
end

TS.s           = s;
TS.component   = component;
TS.intra       = intra;
TS.inter       = inter;
TS.inter_total = sum(inter, 1);
TS.pairs       = pairs;
TS.check       = check;
TS.freqs       = GED.frex(freq_indices);
TS.comp        = comp_idx;
TS.srate       = GED.srate;
end
