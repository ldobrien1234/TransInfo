function [TI_frames_pca,pca_xind,pca_yind]=symm_dynamics()
%This function creates and saves a plot that tracks changes in approximate
%symmetries over the images. The images are ordered alphabetically (ASCII)

    folder = fullfile(pwd, 'data');   % subfolder inside pwd
    files = dir(fullfile(folder,'*_data.mat'));
    files = files(~[files.isdir]); %keep only files, not directories

    %% Iterate through files and plot the minima. Store TI for movie and fPCA
    
    TI_frames=cell(1,numel(files));
    TI_frames_pca=cell(1,numel(files));
    pca_xind=cell(1,numel(files));
    pca_yind=cell(1,numel(files));
    figure;
    for k = 1:numel(files)
        fname = fullfile(folder, files(k).name);
        load(fname, 'TI', 'transRange','scaleRange');
        TI_frames{k}=TI;
        TI_frames_pca{k}=reshape(TI,[1,numel(TI)]);
        %[x_mesh,y_mesh]=meshgrid(transRange,scaleRange);
        [x_mesh,y_mesh]=meshgrid(1:length(transRange),1:length(scaleRange));
        x_mesh=reshape(x_mesh,[1,numel(x_mesh)]);
        y_mesh=reshape(y_mesh,[1,numel(y_mesh)]);
        pca_xind{k}=x_mesh;
        pca_yind{k}=y_mesh;

        %Find approximate symmetries
        %TI=imgaussfilt(TI,1); %smooth out noise to help find local extrema
        [minimaMask,~] = islocalmax2(-TI, 'MinProminence', 0.75);
        [maximaMask,~] = islocalmax2(TI, 'MinProminence', 0.75);
    
        %Extract coordinates
        [minrows, mincols] = find(minimaMask);
        [maxrows,maxcols] = find(maximaMask);
        
        scatter3(transRange(mincols), scaleRange(minrows),k*ones(1,length(minrows)),'r','filled')
        hold on
        scatter3(transRange(maxcols),scaleRange(maxrows),k*ones(1,length(maxrows)),'b','filled')
        hold on
    end
    hold off
    %xlim([0,800])
    %ylim([0.5,1])
    title("Approximate Symmetry Dynamics")
    xlabel('Translation')
    ylabel('Scaling Factor')

    % Setup the VideoWriter object for saving mp4 video files
    v = VideoWriter('mySurfaceAnimation.mp4', 'MPEG-4');
    v.FrameRate = 5; % Frames per second
    open(v);
    
    % Create the figure
    fig = figure;
    ax = axes;
    
    % Loop through plots
    for k = 1:numel(files)

        surf(ax, transRange,scaleRange, TI_frames{k});
        
        zlim([0 2.5]); % Keep axes constant so it doesn't "jump"
        drawnow;        % Force the figure to update
        
        % Capture the frame and write it
        frame = getframe(fig);
        writeVideo(v, frame);
    end
    
    % Close the file
    close(v);

    %% Interactive surface plot
    fig = figure('Name', 'Surface Explorer');
    ax = axes('Parent', fig);
    curr_frame = 1;
    
    % Initial Plot
    plt = surf(ax, TI_frames{curr_frame});
    zlim([min(TI_frames{1}(:)), max(TI_frames{1}(:))]); % Keep scale consistent
    
    % Add a Slider for "Scrubbing"
    uicontrol('Style', 'slider', 'Min', 1, 'Max', numel(files), ...
              'Value', 1, 'SliderStep', [1/numel(files), 0.2*numel(files)], ...
              'Position', [150 20 300 20], ...
              'Callback', @(src, event) update_plot(round(src.Value)));

    function update_plot(frame_idx)
        plt.ZData = TI_frames{frame_idx};
        title(ax, ['Frame: ', num2str(frame_idx)]);
    end


end