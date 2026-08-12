function results = compare_centers_geomean(target, varargin)
%% Compare flower centers: TI rotation vs. TI reflection vs. weighted geometric mean
%TransInfo_images.m (in RotRef) locates two centers for each flower image:
%   1) the center about which rotational TI is minimized ("rotation center")
%   2) the center about which reflectional TI is minimized ("reflection center")
%This script adds a third, much cheaper estimate of the center: a
%function-weighted geometric mean of the pixel coordinates, where each
%pixel is weighted by weightFun(intensity) and only pixels inside the
%flower domain (intensity > Dthresh, matching transinfo's own domain
%definition) contribute. It then plots all three centers on top of each
%image and reports the pairwise distances between them.
%
%   results = compare_centers_geomean()
%   results = compare_centers_geomean(target)
%   results = compare_centers_geomean(target, 'Name', Value, ...)
%
%target - (optional) one of:
%           - a path to a single image file
%           - a path to a folder, searched recursively for .png/.jpg
%             images (optionally limited/subsampled via 'MaxImages')
%           - a cell array (or string array) of image file paths, to
%             process an explicit, hand-picked subset
%          Relative paths are resolved against the repository root (the
%          parent of this file's folder). Default: 'angiosperms'.
%
%Name-Value options:
%  'MaxImages' - cap on the number of images to process when target is a
%                 folder (ignored for a single file or an explicit list).
%                 Default: Inf (process every image found).
%  'Sample'    - how to pick images when MaxImages truncates the folder
%                 listing: 'first' (default) or 'random'.
%  'Seed'      - random seed for 'Sample','random', for reproducibility.
%  'WeightFun' - function handle applied to each pixel's grayscale value
%                 to produce its geometric-mean weight. Default: @(v) v
%                 (i.e. weight by intensity itself).
%
%results - table with one row per image: the three centers, and the
%           pairwise Euclidean distances between them (in pixels). Also
%           written to experiments/results/center_comparison.csv and
%           plotted alongside each image in experiments/plots.
%
%Examples:
%   compare_centers_geomean('angiosperms', 'MaxImages', 5); %first 5 found
%   compare_centers_geomean('angiosperms', 'MaxImages', 5, ...
%       'Sample', 'random', 'Seed', 1); %5 random images
%   compare_centers_geomean({'angiosperms/early_angiosperms/4_nymphaeaceae.png', ...
%       'angiosperms/eudicots/some_flower.png'}); %hand-picked subset

if nargin < 1 || isempty(target)
    target = 'angiosperms';
end

p = inputParser;
addParameter(p, 'MaxImages', Inf, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Sample', 'first', @(x) any(strcmpi(x, {'first','random'})));
addParameter(p, 'Seed', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'WeightFun', @(v) v, @(x) isa(x, 'function_handle'));
parse(p, varargin{:});
maxImages = p.Results.MaxImages;
sampleMode = lower(p.Results.Sample);
seed = p.Results.Seed;
weightFun = p.Results.WeightFun;

Dthresh = 15;
Pmax = 256;

%Resolve paths relative to the repository root (parent of this file's folder)
thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
addpath(fullfile(repoRoot, 'AuxFunctions'));

if iscell(target) || (isstring(target) && ~isscalar(target))
    %Explicit, hand-picked subset of image files. Collect each dir()-style
    %entry into a cell first and vertcat, since dir() structs carry more
    %fields than a hand-rolled struct would, and struct-array assignment
    %requires matching fields throughout.
    entries = cell(numel(target), 1);
    for i = 1:numel(target)
        entries{i} = resolve_image_entry(target{i}, repoRoot);
    end
    files = vertcat(entries{:});
    fullTop = repoRoot;
else
    targetPath = target;
    if ~(isfile(targetPath) || isfolder(targetPath))
        targetPath = fullfile(repoRoot, target);
    end

    if isfile(targetPath)
        files = dir(targetPath);
        fullTop = files(1).folder;
    elseif isfolder(targetPath)
        fullTop = targetPath;
        pngFiles = dir(fullfile(fullTop, '**', '*.png'));
        jpgFiles = dir(fullfile(fullTop, '**', '*.jpg'));
        files = [pngFiles; jpgFiles];

        %Optionally cap/subsample the folder listing
        n = numel(files);
        if isfinite(maxImages) && maxImages < n
            if strcmp(sampleMode, 'random')
                if ~isempty(seed); rng(seed); end
                idx = randperm(n, maxImages);
            else
                idx = 1:maxImages;
            end
            files = files(idx);
        end
    else
        error('compare_centers_geomean:badTarget', 'Could not find image or folder: %s', target);
    end
end

if isempty(files)
    error('compare_centers_geomean:noImages', 'No .png/.jpg images found under: %s', fullTop);
end

plotFolder = fullfile(fileparts(thisFile), 'plots');
resultsFolder = fullfile(fileparts(thisFile), 'results');
if ~exist(plotFolder, 'dir'); mkdir(plotFolder); end
if ~exist(resultsFolder, 'dir'); mkdir(resultsFolder); end

names = strings(numel(files),1);
rxcAll = zeros(numel(files),1); rycAll = zeros(numel(files),1);
bxcAll = zeros(numel(files),1); bycAll = zeros(numel(files),1);
gxcAll = zeros(numel(files),1); gycAll = zeros(numel(files),1);
dRotGeo = zeros(numel(files),1); dRefGeo = zeros(numel(files),1); dRotRef = zeros(numel(files),1);

for k = 1:numel(files)
    imgPath = fullfile(files(k).folder, files(k).name);
    [~, name, ~] = fileparts(files(k).name);
    fprintf('Processing %s ...\n', imgPath);

    imdata0 = imread(imgPath);
    imdata = im2gray(imdata0) + 1;

    [rxc, ryc] = find_rotation_center(imdata, Dthresh, Pmax);
    [bxc, byc] = find_reflection_center(imdata, Dthresh, Pmax);
    [gxc, gyc] = weighted_geomean_center(imdata, Dthresh, weightFun);

    names(k) = name;
    rxcAll(k) = rxc; rycAll(k) = ryc;
    bxcAll(k) = bxc; bycAll(k) = byc;
    gxcAll(k) = gxc; gycAll(k) = gyc;
    dRotGeo(k) = hypot(rxc-gxc, ryc-gyc);
    dRefGeo(k) = hypot(bxc-gxc, byc-gyc);
    dRotRef(k) = hypot(rxc-bxc, ryc-byc);

    %Plot the image with all three centers marked
    h = figure('Visible', 'off');
    ax = axes(h); hold(ax, 'on');
    image(ax, imdata0);
    axis(ax, 'image'); axis(ax, 'off');
    set(ax, 'YDir', 'reverse');
    plot(ax, rxc, ryc, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(ax, bxc, byc, 'c*', 'MarkerSize', 8, 'LineWidth', 1.5);
    plot(ax, gxc, gyc, 'gx', 'MarkerSize', 8, 'LineWidth', 1.5);
    legend(ax, {'rotation center', 'reflection center', 'weighted geometric mean'}, ...
        'Location', 'southoutside', 'TextColor', 'black', 'Color', 'white', 'Box', 'off');
    title(ax, name, 'Interpreter', 'none');
    hold(ax, 'off');

    relativeFolder = erase(files(k).folder, fullTop);
    outFolder = fullfile(plotFolder, relativeFolder);
    if ~exist(outFolder, 'dir'); mkdir(outFolder); end
    exportgraphics(h, fullfile(outFolder, [char(name), '_centers.png']));
    close(h);
end

results = table(names, rxcAll, rycAll, bxcAll, bycAll, gxcAll, gycAll, dRotGeo, dRefGeo, dRotRef, ...
    'VariableNames', {'image', 'rot_x', 'rot_y', 'ref_x', 'ref_y', 'geomean_x', 'geomean_y', ...
    'dist_rot_geomean', 'dist_ref_geomean', 'dist_rot_ref'});

writetable(results, fullfile(resultsFolder, 'center_comparison.csv'));

fprintf('\nMean distance rotation-center to geometric mean: %.2f px\n', mean(dRotGeo, 'omitnan'));
fprintf('Mean distance reflection-center to geometric mean: %.2f px\n', mean(dRefGeo, 'omitnan'));
fprintf('Mean distance rotation-center to reflection-center: %.2f px\n', mean(dRotRef, 'omitnan'));

end


function entry = resolve_image_entry(t, repoRoot)
%Resolve a single target (path relative to repoRoot, or absolute) to a
%dir()-style struct with .folder and .name, for use in an explicit list
%of hand-picked images.
tp = char(t);
if ~isfile(tp)
    tp = fullfile(repoRoot, tp);
end
if ~isfile(tp)
    error('compare_centers_geomean:badTarget', 'Could not find image file: %s', t);
end
d = dir(tp);
entry = d(1);
end


function [rxc, ryc] = find_rotation_center(imdata, Dthresh, Pmax)
%Grid search over candidate centers and rotation angles for the center
%that minimizes TI under rotation. Mirrors the "Find center by rotation
%about x,y" section of TransInfo_images.m.
[M, N] = size(imdata);
pyc = ceil(M/2);
pxc = ceil(N/2);
xcoord = (pxc-floor(M/4)):(pxc+ceil(M/4));
ycoord = (pyc-floor(N/4)):(pyc+ceil(N/4));

NXmax = 21;
NYmax = 21;
xshifts = round(linspace(min(xcoord), max(xcoord), NXmax));
yshifts = round(linspace(min(ycoord), max(ycoord), NYmax));

dTh = pi/6;
thetavals = dTh:dTh:(2*pi-dTh);

%Candidates whose transformed domain overlaps less than COVthresh of the
%original flower's domain are rejected: transinfo's TI is an *average*
%over that overlap, so a wrong-but-far center that rotates most of the
%flower out of its own footprint can leave only a handful of leftover
%pixels whose average divergence is spuriously small, mimicking a true
%symmetry center. See transinfo.m's coverage output for details.
COVthresh = 0.7;

TIcenter = zeros(NXmax, NYmax, length(thetavals));
COVcenter = zeros(NXmax, NYmax, length(thetavals));
for llx = 1:NXmax
    bx = xshifts(llx);
    for kky = 1:NYmax
        by = yshifts(kky);
        Atb = [1,0,0; 0,1,0; -bx,-by,1];
        AtbI = [1,0,0; 0,1,0; bx,by,1];
        for mth = 1:length(thetavals)
            theta = thetavals(mth);
            Arot = [cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
            Atform = Atb*Arot*AtbI;
            [TIcenter(llx,kky,mth), COVcenter(llx,kky,mth)] = transinfo(imdata, Atform, Dthresh, Pmax);
        end
    end
end

TImasked = TIcenter;
TImasked(COVcenter < COVthresh) = Inf;
arr = min(TImasked, [], 3);
minimum = min(arr(:));
if ~isfinite(minimum)
    %No candidate anywhere in the grid kept enough of the flower's domain
    %(e.g. a very small or off-center flower); fall back to the
    %unmasked search rather than erroring out.
    warning('find_rotation_center:noCoverage', ...
        'No candidate center reached %.0f%% domain coverage; falling back to unmasked TI.', COVthresh*100);
    arr = min(TIcenter, [], 3);
    minimum = min(arr(:));
end
[xidmax, yidmax] = find(arr==minimum, 1);
rxc = xshifts(xidmax);
ryc = yshifts(yidmax);
end


function [bxc, byc] = find_reflection_center(imdata, Dthresh, Pmax)
%Search over candidate axis positions for the center that minimizes TI
%under reflection across the vertical/horizontal axis. Mirrors the
%"translate and reflect x,y to find center of TI" section of
%TransInfo_images.m.
[M, N] = size(imdata);
pyc = ceil(M/2);
pxc = ceil(N/2);
xcoord = (pxc-floor(M/4)):(pxc+ceil(M/4));
ycoord = (pyc-floor(N/4)):(pyc+ceil(N/4));

NXmax = 21;
NYmax = 21;
xshifts = round(linspace(min(xcoord), max(xcoord), NXmax));
yshifts = round(linspace(min(ycoord), max(ycoord), NYmax));

Arx = [-1,0,0; 0,1,0; 0,0,1]; %reflect across y-axis (x -> -x)
Ary = [1,0,0; 0,-1,0; 0,0,1]; %reflect across x-axis (y -> -y)

%See find_rotation_center for why low-coverage candidates are rejected.
COVthresh = 0.7;

TIx = zeros(1, NXmax);
COVx = zeros(1, NXmax);
for llx = 1:NXmax
    bx = xshifts(llx);
    Atb = [1,0,0; 0,1,0; -bx,0,1];
    AtbI = [1,0,0; 0,1,0; bx,0,1];
    [TIx(llx), COVx(llx)] = transinfo(imdata, Atb*Arx*AtbI, Dthresh, Pmax);
end

TIy = zeros(1, NYmax);
COVy = zeros(1, NYmax);
for kky = 1:NYmax
    by = yshifts(kky);
    Atb = [1,0,0; 0,1,0; 0,-by,1];
    AtbI = [1,0,0; 0,1,0; 0,by,1];
    [TIy(kky), COVy(kky)] = transinfo(imdata, Atb*Ary*AtbI, Dthresh, Pmax);
end

bxc = xshifts(min_with_endpoints(TIx, COVx, COVthresh));
byc = yshifts(min_with_endpoints(TIy, COVy, COVthresh));
end


function idx = min_with_endpoints(TIvec, COVvec, COVthresh)
%Center will be at a minimum of TIvec. Local minima (peaks of -TIvec) plus
%the two endpoints are the candidates; return the index of the smallest,
%rejecting candidates whose domain coverage falls below COVthresh (see
%find_rotation_center) unless none qualify, in which case fall back to
%the unmasked search.
TImasked = TIvec;
TImasked(COVvec < COVthresh) = Inf;
if ~isfinite(min(TImasked))
    warning('min_with_endpoints:noCoverage', ...
        'No candidate reached %.0f%% domain coverage; falling back to unmasked TI.', COVthresh*100);
    TImasked = TIvec;
end
[pks, locs] = findpeaks(-TImasked);
pks = horzcat(TImasked(1), -pks, TImasked(end));
locs = horzcat(1, locs, length(TImasked));
[~, n] = min(pks);
idx = locs(n);
end


function [gxc, gyc] = weighted_geomean_center(imdata, Dthresh, weightFun)
%Function-weighted geometric mean of the pixel coordinates: each pixel
%(x,y) inside the flower domain (grayscale intensity > Dthresh, matching
%transinfo's own domain definition) is weighted by weightFun(intensity).
%Coordinates are 1-based pixel indices, so they are always positive and
%log() is well-defined.
[M, N] = size(imdata);
[xx, yy] = meshgrid(1:N, 1:M);

intensity = double(imdata);
w = weightFun(intensity);
w(intensity <= Dthresh) = 0;

Wsum = sum(w(:));
if Wsum <= 0
    gxc = NaN;
    gyc = NaN;
    warning('weighted_geomean_center:emptyDomain', ...
        'No pixels above Dthresh (or all weights <= 0); returning NaN.');
    return;
end

gxc = exp(sum(w(:) .* log(xx(:))) / Wsum);
gyc = exp(sum(w(:) .* log(yy(:))) / Wsum);
end
