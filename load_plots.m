%Loads interactive plots
clear

%Adding path for auxiliary functions
addpath(fullfile(pwd,'AuxFunctions'));

%Get full path of selected plot
[filename,plotPath]=uigetfile(pwd,"Please select the .fig plot to display");
%Make plot path relative to the primary working directory
plotPath=erase(plotPath, [pwd,filesep]); 

separatorIndices = strfind(plotPath, filesep);
sepIdx=separatorIndices(1); %Index of first filesep

subFolder=plotPath(1:sepIdx); %includes filesep
plotPath=plotPath(sepIdx+1:end-1); %excludes filesep at beginning and end

%Get top folder
sepFolders=strsplit(plotPath,filesep);
topFolder=sepFolders(1);
topFolder=topFolder{1};

%Cases for the different types of plots that could be selected
if strcmp(topFolder,'pcaPlots') %PCA by class
    dataFile=[fullfile(subFolder,'pcaData',erase(plotPath,[topFolder,filesep])),'.mat'];

    load(dataFile,'filenames');
    fig=openfig(fullfile(subFolder,plotPath,filename),'reuse');

    fig.UserData=filenames;
    
    % Reattach data cursor behavior
    dcm = datacursormode(fig);
    dcm.Enable = 'on';
    dcm.UpdateFcn = @(src, event) showFilename(event, fig);

    % Turn off interpreter for all data tips
    set(findall(fig,'Type','DataTip'),'Interpreter','none');
    
    %Force figure to open
    figure(fig);

elseif strcmp(topFolder,'all_pcaPlots')

    load(fullfile(subFolder,plotPath,'data.mat'),'cumFileNames');
    fig=openfig(fullfile(subFolder,plotPath,filename),'reuse');

    fig.UserData=cumFileNames;
    
    % Reattach data cursor behavior
    dcm = datacursormode(fig);
    dcm.Enable = 'on';
    dcm.UpdateFcn = @(src, event) showFilename(event, fig);

    % Turn off interpreter for all data tips
    set(findall(fig,'Type','DataTip'),'Interpreter','none');
    
    %Force figure to open
    figure(fig);

elseif strcmp(topFolder, 'gwtau_plots')

    load(fullfile(subFolder,plotPath,'data.mat'),'cumFileNames');
    fig=openfig(fullfile(subFolder,plotPath,filename),'reuse');

    fig.UserData=cumFileNames;
    
    % Reattach data cursor behavior
    dcm = datacursormode(fig);
    dcm.Enable = 'on';
    dcm.UpdateFcn = @(src, event) showFilename(event, fig);

    % Turn off interpreter for all data tips
    set(findall(fig,'Type','DataTip'),'Interpreter','none');
    
    %Force figure to open
    figure(fig);


elseif strcmp(topFolder, 'subset_pcaPlots')

    load(fullfile(subFolder,'subset_pcaData','subset_pcaData.mat'),'sorted_labels');
    fig=openfig(fullfile(subFolder,plotPath,filename),'reuse');

    fig.UserData=sorted_labels;
    
    % Reattach data cursor behavior
    dcm = datacursormode(fig);
    dcm.Enable = 'on';
    dcm.UpdateFcn = @(src, event) showFilename(event, fig);

    % Turn off interpreter for all data tips
    set(findall(fig,'Type','DataTip'),'Interpreter','none');
    
    %Force figure to open
    figure(fig);

end
