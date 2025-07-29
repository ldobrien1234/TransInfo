%Loads interactive plots
clear

%Get full path of selected plot
[filename,plotPath]=uigetfile(pwd,"Please select the .fig plot to display");
%Make plot path relative to the primary working directory
plotPath=erase(plotPath, [pwd,filesep]); 


%Get filename then topfolder
currentPath = plotPath;
while true
    [parentDir, folderName, ~] = fileparts(currentPath);
    if isempty(parentDir) || strcmp(parentDir, folderName) % Check for root directory
        topFolder = folderName;
        break;
    end
    topFolder = folderName; % Store the current folder name
    currentPath = parentDir;
end

plotPath=plotPath(1:end-1); %Remove / at end

%Cases for the different types of plots that could be selected
if strcmp(topFolder,'pcaPlots') %PCA by class
    dataFile=[fullfile('pcaData',erase(plotPath,[topFolder,filesep])),'.mat'];

    load(dataFile,'filenames');
    fig=openfig(fullfile(plotPath,filename),'reuse');

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

    load(fullfile(plotPath,'data.mat'),'cumFileNames');
    fig=openfig(fullfile(plotPath,filename),'reuse');

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

    load(fullfile(plotPath,'data.mat'),'cumFileNames');
    fig=openfig(fullfile(plotPath,filename),'reuse');

    fig.UserData=cumFileNames;
    
    % Reattach data cursor behavior
    dcm = datacursormode(fig);
    dcm.Enable = 'on';
    dcm.UpdateFcn = @(src, event) showFilename(event, fig);

    % Turn off interpreter for all data tips
    set(findall(fig,'Type','DataTip'),'Interpreter','none');
    
    %Force figure to open
    figure(fig);

end
