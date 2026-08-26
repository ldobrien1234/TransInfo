
function TI=transinfoUnion(imdata,Atform,Dthresh,Pmax)
%imdata - original image
%Atform - Matrix defining 2D affine transformation
%Dthresh - anything above min(Dthresh,0) of original image will be counted in domain
%Pmax - normalize by TI Pmax. if Pmax=0, normalize by max value of imdata.

tform = affine2d(Atform);

%Keeps output image domain the same as the original image
%Default xlim and ylim for global coordinates are [1.5,ncols+0.5] and
%[1.5,nrows+0.5] (origin at bottom left corner of image)
inDomain=imref2d(size(imdata));
%full option maintains the full bounding box around the transformed image
[~,outDomain] = imwarp(imdata,inDomain,tform);


%Find the intersection of the input and output domains
% Calculate the maximum start and minimum end for both axes
new_XWorldLimits = [min(inDomain.XWorldLimits(1), outDomain.XWorldLimits(1)), ...
                    max(inDomain.XWorldLimits(2), outDomain.XWorldLimits(2))];

new_YWorldLimits = [min(inDomain.YWorldLimits(1), outDomain.YWorldLimits(1)), ...
                    max(inDomain.YWorldLimits(2), outDomain.YWorldLimits(2))];

outputSize=[round(new_XWorldLimits(2)-new_XWorldLimits(1)),...
    round(new_YWorldLimits(2)-new_YWorldLimits(1))];
unionDomain = imref2d(outputSize, new_XWorldLimits, new_YWorldLimits);

%Changing domain of original image and transform Timdata with new bounds
Timdata=imwarp(imdata,inDomain,tform,'OutputView',unionDomain);
imdata=imwarp(imdata,inDomain,affine2d(),'OutputView',unionDomain);

%Set zero values equal to 1
imdata(imdata==0)=1;
Timdata(Timdata==0)=1;

Dim  = double(imdata);
DTim = double(Timdata);

if Pmax==0
    Pmax = max(max(Dim));
end

integrand = Dim.*log(Dim./DTim);

%Ensure no blowups (but shouldn't need)
TI = sum(integrand(isfinite(integrand)));

%normalization
Acount=sum(sum(isfinite(integrand)));
TI=TI/Acount/Pmax;

