function C = hyperness_contributions(GED, X, nchan_per_brain, freq_indices, comp_idx)
% HYPERNESS_CONTRIBUTIONS  Split the power of a multi-brain GED component into
% within-brain (intra) and between-brain (inter) contributions.
%
%   C = hyperness_contributions(GED, X, nchan_per_brain, freq_indices, comp_idx)
%
% INPUT
%   GED             struct with fields .evecs [nChan x nComps x nFreq],
%                   .frex, .srate, .fwidth (see hyperness_pipeline)
%   X               [nChanTotal x nTime] the data the GED was computed on
%   nchan_per_brain e.g. [22 22] for a dyad, [22 22 22] for a triad.
%                   Brains are assumed to be stacked in this order.
%   freq_indices    indices into GED.frex defining the band of interest
%   comp_idx        component to decompose (scalar)
%
% METHOD
%   At each frequency the component power decomposes exactly:
%
%     v' S v = SUM_a v_a' S_aa v_a  +  SUM_{a<b} 2 v_a' S_ab v_b
%              \_____ intra _____/     \________ inter ________/
%
%   Intra terms are non-negative; inter terms are signed (positive = the two
%   brains contribute in phase, negative = anti phase). The decomposition is
%   verified against v'*S*v at every frequency.
%
%   The raw contributions are averaged over the band first and normalised
%   afterwards, i.e. a power-weighted band average (same as the original
%   dyadic implementation).
%
% OUTPUT
%   C.intra          [nBrain x 1]        band-averaged intra ratios
%   C.inter          [nBrain x nBrain]   band-averaged inter ratios; the pair
%                                        value is stored in both (a,b) and
%                                        (b,a), diagonal is zero
%   C.intra_total, C.inter_total         scalars, sum to 1
%   C.intra_raw, C.inter_raw             unnormalised band averages
%   C.intra_ratio_f  [nBrain x nFreq]    per-frequency ratios, to check
%   C.inter_ratio_f  [nBrain x nBrain x nFreq]   consistency inside the band
%   C.total_f        [1 x nFreq]         v'*S*v at each frequency
%   C.inter_chan     {nBrain x nBrain}   inter_chan{a,b} = 2*(v_a v_b') .* S_ab
%                                        normalised and band-averaged, so
%                                        sum(inter_chan{a,b}(:)) == C.inter(a,b)
%   C.freqs, C.comp
%
% DEPENDENCY: filterFGx.m (FREQ-NESS external functions)

nBrain = numel(nchan_per_brain);
nFreq  = numel(freq_indices);

% row indices of each brain inside X
idx = cell(1, nBrain);
o = 0;
for a = 1:nBrain
    idx{a} = (o + 1):(o + nchan_per_brain(a));
    o = o + nchan_per_brain(a);
end
if o ~= size(X, 1)
    error('sum(nchan_per_brain) = %d does not match size(X,1) = %d.', o, size(X, 1));
end

intra_f = zeros(nBrain, nFreq);
inter_f = zeros(nBrain, nBrain, nFreq);
total_f = zeros(1, nFreq);

inter_chan_sum = cell(nBrain, nBrain);
for a = 1:nBrain
    for b = (a+1):nBrain
        inter_chan_sum{a, b} = zeros(nchan_per_brain(a), nchan_per_brain(b));
    end
end

for k = 1:nFreq

    freq_idx = freq_indices(k);

    % ---- narrow-band covariance at this frequency ----
    narrowData = filterFGx(X, GED.srate, GED.frex(freq_idx), GED.fwidth(freq_idx), 0);
    covS = cov(narrowData');

    % ---- eigenvector of the component ----
    v = GED.evecs(:, comp_idx, freq_idx);

    % ---- block decomposition ----
    for a = 1:nBrain
        intra_f(a, k) = real(v(idx{a})' * covS(idx{a}, idx{a}) * v(idx{a}));
    end

    for a = 1:nBrain
        for b = (a+1):nBrain
            % channel-by-channel contribution of this pair
            E = 2 * real((v(idx{a}) * v(idx{b})') .* covS(idx{a}, idx{b}));
            inter_f(a, b, k) = sum(E(:));
            inter_f(b, a, k) = inter_f(a, b, k);
            inter_chan_sum{a, b} = inter_chan_sum{a, b} + E / nFreq;
            % legacy view (spatial filter geometry only, not a contribution):
            %   inter_chan_sum{a,b} = inter_chan_sum{a,b} + (v(idx{a})*v(idx{b})')/nFreq;
        end
    end

    % ---- check the decomposition against the full quadratic form ----
    total_f(k)  = sum(intra_f(:, k)) + sum(sum(triu(inter_f(:, :, k), 1)));
    full_contrib = real(v' * covS * v);
    if abs(total_f(k) - full_contrib) > 1e-10 * max(1, abs(full_contrib))
        warning('hyperness:decomposition', ...
            'decomposition does not match the full quadratic form at %.4f Hz: difference = %.2e', ...
            GED.frex(freq_idx), abs(total_f(k) - full_contrib));
    end
end

% ---- band average: average the raw contributions, then normalise ----
intra_raw = mean(intra_f, 2);
inter_raw = mean(inter_f, 3);
den = sum(intra_raw) + sum(sum(triu(inter_raw, 1)));

C.intra       = intra_raw / den;
C.inter       = inter_raw / den;
C.intra_total = sum(C.intra);
C.inter_total = sum(sum(triu(C.inter, 1)));
C.intra_raw   = intra_raw;
C.inter_raw   = inter_raw;

% per-frequency ratios, to check whether the band is homogeneous
C.intra_ratio_f = intra_f ./ repmat(total_f, nBrain, 1);
C.inter_ratio_f = inter_f ./ repmat(reshape(total_f, 1, 1, nFreq), nBrain, nBrain, 1);
C.total_f       = total_f;

C.inter_chan = cell(nBrain, nBrain);
for a = 1:nBrain
    for b = (a+1):nBrain
        C.inter_chan{a, b} = inter_chan_sum{a, b} / den;
    end
end

C.freqs = GED.frex(freq_indices);
C.comp  = comp_idx;
end
