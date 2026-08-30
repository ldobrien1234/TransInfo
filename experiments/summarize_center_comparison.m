function summarize_center_comparison(target)
%% Summarize/plot the pairwise-distance columns of an existing center_comparison.csv
%compare_centers_centroid.m and compare_centers_centroid_union.m each write
%a center_comparison.csv (to resultsAll/ and resultsUnion/ respectively);
%this script derives the rest of that folder's output from it: a
%per-metric summary-statistics CSV and one histogram per metric.
%
%Only the percentage-of-image-diagonal columns (*_pct) are summarized and
%plotted, not the raw pixel columns: raw pixel distances aren't comparable
%across differently-sized images, so the percentage columns are the
%meaningful ones. The three histograms share one x-axis range and one
%y-axis range (rather than each autoscaling independently, as the
%previous histograms did), so their shapes are directly comparable.
%
%   summarize_center_comparison()
%   summarize_center_comparison(folder)
%   summarize_center_comparison({folder1, folder2, ...})
%
%target - (optional) one of:
%           - a path to a folder containing center_comparison.csv (e.g.
%             'resultsAll' or 'resultsUnion'), resolved against this
%             file's folder if not found as given
%           - a cell array of such folders
%          Default: {'resultsAll', 'resultsUnion'}.
%
%For each folder, writes:
%   <folder>/center_comparison_pct_summary.csv - one row per *_pct metric:
%       n, mean, median, standard deviation, 90th/95th percentile, and the
%       percentage of images under 5%/10%.
%   <folder>/dist_rot_centroid_pct_hist.png
%   <folder>/dist_ref_centroid_pct_hist.png
%   <folder>/dist_rot_ref_pct_hist.png
%and deletes any stale raw-pixel *_hist.png files (dist_rot_centroid_hist.png
%etc.) left over from an older run, since this script no longer produces them.
%
%Examples:
%   summarize_center_comparison();              % both resultsAll and resultsUnion
%   summarize_center_comparison('resultsAll');   % just one folder

if nargin < 1 || isempty(target)
    target = {'resultsAll', 'resultsUnion'};
end
if ischar(target) || (isstring(target) && isscalar(target))
    target = {char(target)};
end

thisFolder = fileparts(mfilename('fullpath'));

for i = 1:numel(target)
    folder = char(target{i});
    if ~isfolder(folder)
        folder = fullfile(thisFolder, folder);
    end
    summarize_one_folder(folder);
end

end


function summarize_one_folder(folder)
csvPath = fullfile(folder, 'center_comparison.csv');
if ~isfile(csvPath)
    error('summarize_center_comparison:missingCsv', ...
        'No center_comparison.csv found in: %s', folder);
end
fprintf('Summarizing %s ...\n', csvPath);
results = readtable(csvPath);

distFields = {'dist_rot_centroid_pct', 'dist_ref_centroid_pct', 'dist_rot_ref_pct'};
distLabels = {'rotation center to weighted centroid', ...
    'reflection center to weighted centroid', 'rotation center to reflection center'};
nMetrics = numel(distFields);

%% Summary statistics
metric = strings(nMetrics,1);
pair = strings(nMetrics,1);
n = zeros(nMetrics,1);
mean_pct = zeros(nMetrics,1);
median_pct = zeros(nMetrics,1);
sd_pct = zeros(nMetrics,1);
p90_pct = zeros(nMetrics,1);
p95_pct = zeros(nMetrics,1);
pct_images_under_5pct = zeros(nMetrics,1);
pct_images_under_10pct = zeros(nMetrics,1);

vals = cell(nMetrics,1);
for i = 1:nMetrics
    v = results.(distFields{i});
    v = v(~isnan(v));
    vals{i} = v;

    metric(i) = distFields{i};
    pair(i) = distLabels{i};
    n(i) = numel(v);
    mean_pct(i) = round(mean(v), 2);
    median_pct(i) = round(median(v), 2);
    sd_pct(i) = round(std(v), 2);
    p90_pct(i) = round(percentile_linear(v, 90), 2);
    p95_pct(i) = round(percentile_linear(v, 95), 2);
    pct_images_under_5pct(i) = round(100 * sum(v < 5) / numel(v), 1);
    pct_images_under_10pct(i) = round(100 * sum(v < 10) / numel(v), 1);
end

summary = table(metric, pair, n, mean_pct, median_pct, sd_pct, p90_pct, p95_pct, ...
    pct_images_under_5pct, pct_images_under_10pct);
writetable(summary, fullfile(folder, 'center_comparison_pct_summary.csv'));

%% Histograms: one shared x-axis and y-axis range across all three metrics
%so their shapes are directly comparable, rather than each autoscaling
%independently.
allVals = vertcat(vals{:});
edges = linspace(0, max(allVals), 31); %30 bins, shared across all three

counts = cell(nMetrics,1);
maxCount = 0;
for i = 1:nMetrics
    counts{i} = histcounts(vals{i}, edges);
    maxCount = max(maxCount, max(counts{i}));
end

for i = 1:nMetrics
    h = figure('Visible', 'off');
    histogram('BinEdges', edges, 'BinCounts', counts{i});
    xlim([edges(1), edges(end)]);
    ylim([0, maxCount * 1.05]);
    xlabel('distance (% of image diagonal)'); ylabel('count');
    title(distLabels{i}, 'Interpreter', 'none');
    exportgraphics(h, fullfile(folder, [distFields{i}, '_hist.png']));
    close(h);
end

%% Clean up stale raw-pixel histograms this script no longer produces
staleFiles = {'dist_rot_centroid_hist.png', 'dist_ref_centroid_hist.png', 'dist_rot_ref_hist.png'};
for i = 1:numel(staleFiles)
    f = fullfile(folder, staleFiles{i});
    if isfile(f); delete(f); end
end

end


function p = percentile_linear(v, pct)
%Linear-interpolation percentile between order statistics (numpy/R-7
%convention): matches the values already recorded in
%center_comparison_pct_summary.csv, which MATLAB's own prctile (an R-5-style
%method) does not reproduce.
v = sort(v(:));
n = numel(v);
if n == 0
    p = NaN;
    return;
end
idx = (pct/100) * (n-1) + 1; %1-based fractional index
lo = floor(idx); hi = ceil(idx);
if lo == hi
    p = v(lo);
else
    frac = idx - lo;
    p = v(lo) + frac * (v(hi) - v(lo));
end
end
