function [y, fwdMean, bwdMean] = ckfilter(x, m, p)
% CKFILTER first order Chung–Kennedy nonlinear forward/backward filter
% (i.e. using one forward and one backward window)
%
%   y = ckfilter(x, m) filters column- or row-vector x with window length m.
%   [y, fwdMean, bwdMean] also returns the forward (future) and backward (past)
%   running means used internally. p is a weighting factor (See Chung &
%   Kennedy (1991))
%
%   y = ckfilter(x, m, 'k', 1.5, 'Pad', 'reflect') lets you tune behavior:
%       'k'    : edge sensitivity (default 1.5). Larger => more averaging,
%                smaller => more aggressive edge preservation.
%       'Pad'  : padding at the ends: 'reflect' (default), 'replicate', or 'none'.
%
% THEORY (practical version)
%   For each index i:
%     - bwdMean(i) = mean(x[i-m+1 : i])           (past window)
%     - fwdMean(i) = mean(x[i : i+m-1])           (future window)
%     - bwdVar_m(i), fwdVar_m(i) = variance of the sample means
%       estimated from running variances divided by m.
%   Let Delta = |fwdMean - bwdMean|.
%   Let scale = sqrt(fwdVar_m + bwdVar_m).
%   If Delta <= k * scale (locally “flat”): do precision-weighted average
%       y = (fwdMean/fwdVar_m + bwdMean/bwdVar_m) / (1/fwdVar_m + 1/bwdVar_m)
%   Else (near an edge): choose the side better explaining x(i):
%       y = argmin_{mu ∈ {fwdMean, bwdMean}} (x(i) - mu)^2 / var_mu
%
% NOTES
%   - m must be a positive integer >= 2.
%   - This is faithful to the CK spirit (forward/backward, adaptive, edge-preserving)
%     and works very well in practice for step-like signals in Gaussian noise.
%   - For matrices, call per column (e.g., y = ckfilter(x(:,j), m)).

%% Validate input vector
if ~isvector(x)
    error('x must be a vector.');
end
x = x(:);                        % work in column form
n = numel(x);
if ~(isscalar(m) && m==round(m) && m>=2)
    error('m (window length) must be an integer >= 2.');
end

%% Pad signal by reflecting on the boundary
pre  = flipud(x(2:min(n, m)));
post = flipud(x(max(1, n-m+1):n-1));
xp   = [pre; x; post];
off  = numel(pre);

%% Running forward and backward means, and sum of squared difference
% Helper: running mean via causal box filter
avgBox = ones(m,1)/m;
box=ones(m,1);
padding=10^(-3); %to avoid division by zero in weighting scheme

% Backward (past) mean at i: mean of last m samples up to i
bwdMean_p = filter(avgBox, 1, xp);

% Forward (future) mean at i: do causal mean on flipped signal, then flip back
fwdMean_p = flipud(filter(avgBox, 1, flipud(xp)));

%Summing over future time points for backwards weight
bwdWeight_p=flipud(filter(box,1,flipud((xp-bwdMean_p).^2))+padding).^(-p);

%Summing over previous time points for forward weight
fwdWeight_p=(filter(box,1,(xp-fwdMean_p).^2)+padding).^(-p);


% Trim padding back to original length
if off>0
    idx = (off+1):(off+n);
    bwdMean = bwdMean_p(idx);
    fwdMean = fwdMean_p(idx);
    bwdWeight = bwdWeight_p(idx);
    fwdWeight = fwdWeight_p(idx);
    %xw      = xp(idx);
else
    bwdMean = bwdMean_p;
    fwdMean = fwdMean_p;
    bwdWeight = bwdWeight_p;
    fwdWeight = fwdWeight_p;
    %xw      = xp;
end

%Normalize weights
total=fwdWeight+bwdWeight;
fwdWeight_n=fwdWeight./total;
bwdWeight_n=bwdWeight./total;

y=fwdWeight_n.*fwdMean+bwdWeight_n.*bwdMean;

%Restore original shape
if isrow(x), y = y.'; fwdMean = fwdMean.'; bwdMean = bwdMean.'; end
end

