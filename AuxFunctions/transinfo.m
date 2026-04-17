function TI=transinfo(imdata,Atform,Dthresh,Pmax,option)
%imdata - original image
%Atform - Matrix defining 2D affine transformation
%Dthresh - domain cutoff: pixels whose value exceeds Dthresh count as
%          flower (i.e. inside the domain D)
%Pmax - normalize by TI Pmax. if Pmax=0, normalize by max value of imdata.
%option - (optional) selects the domain-masking method. Default: 'hard'.
%   'hard' : hard indicator mask, domain = {Dim > Dthresh & DTim > Dthresh}.
%            Fast. Under noise, however, individual pixels cross the
%            threshold one at a time, and at rational rotation angles
%            (pi/2, pi/3, pi, 2pi/3, ...) bilinear interpolation gives many
%            pixels identical interpolation weights, so entire cohorts of
%            pixels cross the threshold simultaneously. This manifests as
%            visible spikes in the TI curve at those angles.
%   'soft' : cubic-Hermite smoothstep weight that ramps from 0 at
%            Dthresh-delta to 1 at Dthresh+delta. The weight is a
%            continuous function of each pixel's value, which in turn is
%            continuous in theta (bilinear interpolation is continuous),
%            so TI is continuous in theta and the rational-angle spikes
%            go away. Default delta = Dthresh/3.
%   struct('method','soft','delta',X) : soft mask with a custom
%            transition half-width X.

if nargin < 5 || isempty(option)
    method = 'hard';
    delta  = max(Dthresh,1)/3;
elseif ischar(option) || isstring(option)
    method = char(option);
    delta  = max(Dthresh,1)/3;
elseif isstruct(option) && isfield(option,'method')
    method = option.method;
    if isfield(option,'delta') && ~isempty(option.delta)
        delta = option.delta;
    else
        delta = max(Dthresh,1)/3;
    end
else
    error('transinfo: unrecognized option argument (use ''hard'', ''soft'', or a struct).');
end

[M,N] = size(imdata);
tform = affine2d(Atform);

%Keeps output image domain the same as the original image
%Default xlim and ylim for global coordinates are [1.5,ncols+0.5] and
%[1.5,nrows+0.5] (origin at bottom left corner of image)
inDomain=imref2d(size(imdata));
%full option maintains the full bounding box around the transformed image
[Timdata,outDomain] = imwarp(imdata,inDomain,tform,'OutputView',inDomain);


%Find the intersection of the input and output domains
% Calculate the maximum start and minimum end for both axes
new_XWorldLimits = [max(inDomain.XWorldLimits(1), outDomain.XWorldLimits(1)), ...
                    min(inDomain.XWorldLimits(2), outDomain.XWorldLimits(2))];

new_YWorldLimits = [max(inDomain.YWorldLimits(1), outDomain.YWorldLimits(1)), ...
                    min(inDomain.YWorldLimits(2), outDomain.YWorldLimits(2))];

intersection = imref2d([M,N], new_XWorldLimits, new_YWorldLimits);

Timdata=imwarp(Timdata,outDomain,affine2d(),'OutputView',intersection);
%By default, pixels that in the transformed image outside of interpolation
%are set to zero. These are ignored in the integral since they are outside
%of Dtilde

Dim  = double(imdata);
DTim = double(Timdata);

if Pmax==0
    Pmax = max(max(Dim));
end

threshold = max(Dthresh,0);

switch lower(method)
    case 'hard'
        %Hard-indicator domain. See header comment for why this produces
        %rational-angle spikes under noise.
        domain = (Dim > threshold) & (DTim > threshold);
        Acount = sum(domain(:));

        if Acount == 0
            TI = 0;
            return;
        end

        integrand = Dim .* log(Dim ./ DTim);
        TI = sum(integrand(domain));
        TI = TI / Acount / Pmax;

    case 'soft'
        %Soft cubic-Hermite smoothstep weight. Requires threshold-delta>=0
        %so that pixels with DTim=0 (zero-padded boundary) receive exactly
        %zero weight and the (finite-value * zero-weight) product at the
        %boundary is well-defined. Clamp delta so this holds.
        d = min(delta, threshold);

        wD  = soft_weight(Dim,  threshold, d);
        wDT = soft_weight(DTim, threshold, d);
        w   = wD .* wDT;
        Wsum = sum(w(:));

        if Wsum == 0
            TI = 0;
            return;
        end

        %Only evaluate log where DTim>0; weight at DTim=0 is zero by
        %construction (since threshold-d >= 0 and the smoothstep vanishes
        %there), so the unused entries contribute nothing.
        integrand    = zeros(size(Dim));
        safe         = DTim > 0;
        integrand(safe) = Dim(safe) .* log(Dim(safe) ./ DTim(safe));

        TI = sum(w(:) .* integrand(:)) / Wsum / Pmax;

    otherwise
        error('transinfo: unknown method ''%s'' (use ''hard'' or ''soft'').', method);
end

end


function w = soft_weight(x, threshold, delta)
%Cubic Hermite smoothstep:
%   w = 0           for x <= threshold - delta
%   w = 1           for x >= threshold + delta
%   w = smooth      for threshold - delta < x < threshold + delta
%with w'=0 at both endpoints so the composition w(Dim(theta)) is C^1 in
%theta away from a zero-measure set.
if delta <= 0
    %Degenerate to hard indicator
    w = double(x > threshold);
    return;
end
t = (x - (threshold - delta)) / (2*delta);
t = max(0, min(1, t));
w = t .* t .* (3 - 2*t);
end
