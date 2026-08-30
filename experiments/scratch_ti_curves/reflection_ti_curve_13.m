%% Reflection TI curves for angiosperms/early_angiosperms/13.png
%For flower image 13, computes TI as a function of the reflection axis
%angle theta (the flower is reflected about a line through a fixed
%center, at angle theta to the x-axis; theta swept 0:360 degrees),
%about two different centers:
%   1) the reflection center (bxc,byc) = (247,437), i.e. the center that
%      minimizes TI under reflection about the horizontal/vertical axes,
%      as found by find_reflection_center in compare_centers_centroid.m
%      and marked with a cyan star in plotsAll/early_angiosperms/13_centers.png
%      (and recorded in resultsAll/center_comparison.csv, row image=13).
%   2) the weighted centroid (cxc,cyc), recomputed here the same way as
%      weighted_centroid_center in compare_centers_centroid.m.
%TI is computed both ways (transinfo.m, which normalizes by the
%intersection of the original and transformed domains) about both
%centers, and additionally with transinfoUnion.m (which normalizes by the
%union instead) about the reflection center, to compare the two
%normalizations.
%
%Reflection-about-angle-theta transform mirrors the "rotate and reflect
%about center" section of RotRef/TransInfo_images.m: bring the axis at
%angle theta to the x-axis (rotate by -theta), reflect y -> -y, then
%rotate back by theta -- all about the given center.
%
%Saves scratch_ti_curves/13_reflection_ti_curve.mat (thetavals, both
%centers, and the three TI curves) and a plot,
%scratch_ti_curves/13_reflection_ti_curve.png.

%Resolve paths relative to the repository root (parent of experiments/,
%which is this script's grandparent folder)
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(fileparts(thisFile)));
addpath(fullfile(repoRoot, 'AuxFunctions'));

imgPath = fullfile(repoRoot, 'angiosperms', 'early_angiosperms', '13.png');
imdata0 = imread(imgPath);
imdata = im2gray(imdata0) + 1;

Dthresh = 15;
Pmax = 256;

%Reflection center, as found by find_reflection_center in
%compare_centers_centroid.m and recorded in
%resultsAll/center_comparison.csv (row image=13; also the cyan star in
%plotsAll/early_angiosperms/13_centers.png).
bxc = 247;
byc = 437;

%Weighted centroid: sum(w.*coord)/sum(w) over pixels inside the flower
%domain (intensity > Dthresh), weighted by intensity itself -- same
%definition as weighted_centroid_center in compare_centers_centroid.m.
[cxc, cyc] = weighted_centroid_center(imdata, Dthresh, @(v) v);

Nthmax = 361;
thetavals = linspace(0, 360, Nthmax) * pi/180;

Ary = [1,0,0; 0,-1,0; 0,0,1]; %reflect across x-axis (y -> -y)

TI_ref = zeros(Nthmax, 1);       %about reflection center, transinfo
TI_cen = zeros(Nthmax, 1);       %about weighted centroid, transinfo
TI_ref_union = zeros(Nthmax, 1); %about reflection center, transinfoUnion

b_ref = [bxc, byc];
Atb_ref = [1,0,0; 0,1,0; -b_ref(1),-b_ref(2),1];
AtbI_ref = [1,0,0; 0,1,0; b_ref(1),b_ref(2),1];

b_cen = [cxc, cyc];
Atb_cen = [1,0,0; 0,1,0; -b_cen(1),-b_cen(2),1];
AtbI_cen = [1,0,0; 0,1,0; b_cen(1),b_cen(2),1];

for llth = 1:Nthmax
    theta = thetavals(llth);
    Arot = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
    ArotI = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];

    %reflect about the axis at angle theta through the reflection center:
    %bring axis to x-axis, reflect y->-y, rotate back
    tformRef = Atb_ref*ArotI*Ary*Arot*AtbI_ref;
    TI_ref(llth) = transinfo(imdata, tformRef, Dthresh, Pmax);
    TI_ref_union(llth) = transinfoUnion(imdata, tformRef, Dthresh, Pmax);

    %same, but about the weighted centroid
    tformCen = Atb_cen*ArotI*Ary*Arot*AtbI_cen;
    TI_cen(llth) = transinfo(imdata, tformCen, Dthresh, Pmax);
end

outFolder = fileparts(thisFile);
save(fullfile(outFolder, '13_reflection_ti_curve.mat'), ...
    'thetavals', 'bxc', 'byc', 'cxc', 'cyc', 'TI_ref', 'TI_cen', 'TI_ref_union');

h = figure('Visible', 'off');

subplot(2,1,1);
plot(thetavals*180/pi, TI_ref, 'c', 'LineWidth', 1.5); hold on;
plot(thetavals*180/pi, TI_cen, 'g', 'LineWidth', 1.5); hold off;
xlabel('reflection axis angle (degrees)'); ylabel('TI (transinfo)');
legend({'about reflection center', 'about weighted centroid'}, 'Location', 'best');
title('Flower 13: TI under reflection (transinfo)');

subplot(2,1,2);
plot(thetavals*180/pi, TI_ref_union, 'b', 'LineWidth', 1.5);
xlabel('reflection axis angle (degrees)'); ylabel('TI (transinfoUnion)');
legend({'about reflection center'}, 'Location', 'best');
title('Flower 13: TI under reflection (transinfoUnion)');

exportgraphics(h, fullfile(outFolder, '13_reflection_ti_curve.png'));
close(h);

[minRef, iRef] = min(TI_ref);
[minCen, iCen] = min(TI_cen);
[minRefU, iRefU] = min(TI_ref_union);
fprintf('Reflection center: (%.4f, %.4f)\n', bxc, byc);
fprintf('Weighted centroid: (%.4f, %.4f)\n', cxc, cyc);
fprintf('Min TI_ref (transinfo): %.4f at theta=%.1f deg\n', minRef, thetavals(iRef)*180/pi);
fprintf('Min TI_cen (transinfo): %.4f at theta=%.1f deg\n', minCen, thetavals(iCen)*180/pi);
fprintf('Min TI_ref_union (transinfoUnion): %.4f at theta=%.1f deg\n', minRefU, thetavals(iRefU)*180/pi);


function [cxc, cyc] = weighted_centroid_center(imdata, Dthresh, weightFun)
%Function-weighted centroid (weighted arithmetic mean) of the pixel
%coordinates: each pixel (x,y) inside the flower domain (grayscale
%intensity > Dthresh, matching transinfo's own domain definition) is
%weighted by weightFun(intensity), and the center is
%sum(w.*coord)/sum(w) over each coordinate. Mirrors
%weighted_centroid_center in compare_centers_centroid.m.
[M, N] = size(imdata);
[xx, yy] = meshgrid(1:N, 1:M);

intensity = double(imdata);
w = weightFun(intensity);
w(intensity <= Dthresh) = 0;

Wsum = sum(w(:));
if Wsum <= 0
    cxc = NaN;
    cyc = NaN;
    warning('weighted_centroid_center:emptyDomain', ...
        'No pixels above Dthresh (or all weights <= 0); returning NaN.');
    return;
end

cxc = sum(w(:) .* xx(:)) / Wsum;
cyc = sum(w(:) .* yy(:)) / Wsum;
end
