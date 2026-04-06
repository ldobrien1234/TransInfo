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

Acount    = sum(sum(Dim>0 & DTim>0));
integrand = Dim.*log(Dim./DTim);

%Values outside of \tilde{D} are zero for the transformed image;
%they give inf and hence are ignored
TI = sum(integrand(~isinf(integrand)));


TI=TI/Acount/Pmax; %normalization

