function TI=transinfo(imdata,Atform,Dthresh,Pmax,option)
%imdata - original image
%Atform - Matrix defining 2D affine transformation
%Dthresh - anything above min(Dthresh,0) of original image will be counted in domain
%Pmax - normalize by TI Pmax. if Pmax=0, normalize by max value of imdata.

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

%Domain: pixels where both the original and transformed images exceed
%Dthresh. Using Dthresh for DTim (instead of 0) is critical: at angles
%that are exact multiples of pi/2 the rotation aligns with the pixel
%grid so imwarp places boundary pixels at exactly 0 (excluded), whereas
%at nearby angles bilinear interpolation gives tiny positive values like
%0.01 that are included and make Dim*log(Dim/DTim) huge, creating
%artificial jumps in the TI curve.
threshold = max(Dthresh, 0);
domain    = (Dim > threshold) & (DTim > threshold);
Acount    = sum(domain(:));

if Acount == 0
    TI = 0;
    return;
end

integrand = Dim .* log(Dim ./ DTim);
TI = sum(integrand(domain));

TI = TI / Acount / Pmax; %normalization

