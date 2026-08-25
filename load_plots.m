%LOAD_PLOTS Open a saved .fig scatter plot with interactive data tips.
%
%   Run from anywhere; paths resolve relative to this script's own
%   location (TransInfo), so it can load plots from RotRef/ or any other
%   subproject folder alongside it. Prompts you to pick a .fig file, then
%   reopens it with data-cursor mode enabled so that clicking on any
%   point shows the name of the flower/file that point corresponds to.
%
%   Label data is located one of two ways:
%     1. Preferred: a "data.mat" file saved alongside the .fig files
%        (e.g. by ti_pca_all.m), containing a cell array of labels
%        (cumFileNames, filenames, or sorted_labels) in the same order
%        the points were plotted.
%     2. Legacy layouts, where label data lives elsewhere relative to the
%        subproject root that contains the .fig's top-level plot folder:
%           <subproject>/pcaPlots/...        -> <subproject>/pcaData/....mat  ('filenames')
%           <subproject>/subset_pcaPlots/... -> <subproject>/subset_pcaData/subset_pcaData.mat  ('sorted_labels')

clear

%Resolve paths relative to this script's own location, not MATLAB's
%current working directory, so it can be run from anywhere.
scriptDir=fileparts(mfilename('fullpath'));

%Adding path for auxiliary functions
addpath(fullfile(scriptDir,'AuxFunctions'));

%Prompt user to select a .fig file to display
[filename,plotPath]=uigetfile(fullfile(scriptDir,'*.fig'),"Please select the .fig plot to display");

if isequal(filename,0)
    return; %user cancelled the dialog
end

labelVarNames={'cumFileNames','filenames','sorted_labels'};
labels=[];

% --- Preferred: a data.mat co-located with the .fig ---
dataFile=fullfile(plotPath,'data.mat');

if isfile(dataFile)
    dataVars=who('-file',dataFile);
    labelVar=labelVarNames(ismember(labelVarNames,dataVars));
    if ~isempty(labelVar)
        S=load(dataFile,labelVar{1});
        labels=S.(labelVar{1});
    end
end

% --- Legacy layouts, keyed off the subproject/top-level plot folder ---
if isempty(labels)
    relPath=erase(plotPath,[scriptDir,filesep]);
    separatorIndices=strfind(relPath,filesep);

    if ~isempty(separatorIndices)
        sepIdx=separatorIndices(1); %index of first filesep

        subFolder=relPath(1:sepIdx); %e.g. 'RotRef/' -- includes filesep
        subPath=relPath(sepIdx+1:end-1); %e.g. 'pcaPlots/eudicots/rosids' -- excludes leading/trailing filesep

        sepFolders=strsplit(subPath,filesep);
        topFolder=sepFolders{1};

        subFolder=fullfile(scriptDir,subFolder);

        if strcmp(topFolder,'pcaPlots')
            legacyDataFile=[fullfile(subFolder,'pcaData',erase(subPath,[topFolder,filesep])),'.mat'];
            if isfile(legacyDataFile)
                S=load(legacyDataFile,'filenames');
                labels=S.filenames;
            end
        elseif strcmp(topFolder,'subset_pcaPlots')
            legacyDataFile=fullfile(subFolder,'subset_pcaData','subset_pcaData.mat');
            if isfile(legacyDataFile)
                S=load(legacyDataFile,'sorted_labels');
                labels=S.sorted_labels;
            end
        end
    end
end

if isempty(labels)
    error('load_plots:noLabels', ...
        'Could not locate a recognized label variable for %s.',fullfile(plotPath,filename));
end

fig=openfig(fullfile(plotPath,filename),'reuse');
fig.UserData=labels;

% Reattach data cursor behavior
dcm = datacursormode(fig);
dcm.Enable = 'on';
dcm.UpdateFcn = @(src, event) showFilename(event, fig);

% Turn off interpreter for all data tips
set(findall(fig,'Type','DataTip'),'Interpreter','none');

%Force figure to open
figure(fig);
