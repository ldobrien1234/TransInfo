function ti_pca_all()
    
    %Path for functional PCA package
    pacePath="/Users/liamobrien/Documents/PACE_matlab-master/release2.17/PACE";
    addpath(genpath(pacePath));
    %Option to get 2 principal components each time
    opt = setOptions();
    opt.selection_k = 'FVE';
    opt.maxk = 2;
    opt.FVE_threshold = 0.99;

    
    
    %Loop through pca Data file
    pcaDataFolders=genpath("pcaData");
    pcaFoldersList=strsplit(pcaDataFolders,pathsep);


    %initializing
    class_count=0; %initializing "vector" (end+1) indexing still works
    tirotr_all_data=[];
    tirefr_all_data=[];
    tirotb_all_data=[];
    tirefb_all_data=[];

    classNames={}; %To save foldernames for legend
    cumFileNames={}; %Saving filenames for interactive plots

    %Iterating through folder paths
    for i=1:length(pcaFoldersList)
        pcaFolderPath=pcaFoldersList{i};

        if isempty(pcaFolderPath)
            continue;
        end
        
        %lists folder contents with .mat ending
        pca_files=dir(fullfile(pcaFolderPath, '*.mat'));

        if isempty(pca_files)
            continue;
        end
        
        %Number of images in previous classes and current (in loop)
        for i=1:length(pca_files)
            file=pca_files(i);
            load(fullfile(pcaFolderPath,file.name),'tirotr_all','tirefr_all','tirotb_all','tirefb_all','filenames');

            %Rows are TI data for a given flower, columns are theta mesh points
            [rws,~]=size(tirotr_all); %Size of TI data function
            classNames{end+1}=erase(file.name,".mat");
            cumFileNames=[cumFileNames,filenames];
            class_count(end+1)=rws; %Number of flowers in 1st class (early angiosperms)
            tirotr_all_data(end+1:end+rws,:)=tirotr_all;
            tirefr_all_data(end+1:end+rws,:)=tirefr_all;
            tirotb_all_data(end+1:end+rws,:)=tirotb_all;
            tirefb_all_data(end+1:end+rws,:)=tirefb_all;
            %summing previous classes with current for next loop

        end %end loop through files in a folder

    end %end loop through folders in pcaData

    %Converting matrices to cell arrays for functional pca
    tirotr_all_cell=num2cell(tirotr_all_data,2)';
    tirefr_all_cell=num2cell(tirefr_all_data,2)';
    tirotb_all_cell=num2cell(tirotb_all_data,2)';
    tirefb_all_cell=num2cell(tirefb_all_data,2)';

    domainPoints=linspace(0,2*pi,length(tirotr_all_data(1,:)));
    %convert to cell that matches the number of TI curves
    domainPoints=repmat({domainPoints}, numel(tirotr_all_cell), 1)';

    
    class_num=length(class_count)-1; %includes empty class at beginning
    sz=50; %size of points in scatter plots
    colors=lines(class_num); %colors for each class (not including empty class)
    
    tirotr_all_fpca=FPCA(tirotr_all_cell,domainPoints,opt);
    score=getVal(tirotr_all_fpca,'xi_est'); %scores
    fpcs=getVal(tirotr_all_fpca,'phi'); %fractional principal components
    fve=getVal(tirotr_all_fpca,'FVE'); %fraction of variance explained by each pc
    pve=fve(2)*100;
    hrotr=figure('Visible','off'); 
    for class=1:class_num
        %rows in class
        class_rws=sum(class_count(1:class))+1:sum(class_count(1:class+1));
        %Getting first and second principle component for each
        %Coloring different classes different colors
        scatter(score(class_rws,1),score(class_rws,2),sz,colors(class,:),'filled'); 
        hold on
    end
    hold off
    title(append('TI by rotation, Rot ctr.',newline,string(pve),'% variance explained.') );
    legend(classNames,'Interpreter','none'); %add legend for each color
    % print(hrotr,'-djpeg','tirotr_all_angiosperms_pca.jpg')
    h_tirotr_fpc1=figure('Visible','off');
    scatter(domainPoints{1},fpcs(:,1),sz,'filled');
    title("tirotr fpc1")
    h_tirotr_fpc2=figure('Visible','off');
    scatter(domainPoints{1},fpcs(:,2),sz,'filled');
    title("tirotr fpc2")






    tirefr_all_fpca=FPCA(tirefr_all_cell,domainPoints,opt);
    score=getVal(tirefr_all_fpca,'xi_est');
    fpcs=getVal(tirefr_all_fpca,'phi');
    fve=getVal(tirefr_all_fpca,'FVE'); %fraction of variance explained by each pc
    pve=fve(2)*100;
    hrefr=figure('Visible','off');
    for class=1:class_num
        %rows in class
        class_rws=sum(class_count(1:class))+1:sum(class_count(1:class+1));
        %Getting first and second principle component for each
        %Coloring different classes different colors
        scatter(score(class_rws,1),score(class_rws,2),sz,colors(class,:),'filled');
        hold on
    end
    hold off
    title(append('TI by reflection Rot ctr.',newline,string(pve),'% variance explained.'));
    legend(classNames,'Interpreter','none'); %add legend for each color
    % print(hrefr,'-djpeg','tirefr_all_angiosperms_pca.jpg')
    h_tirefr_fpc1=figure('Visible','off');
    scatter(domainPoints{1},fpcs(:,1),sz,'filled');
    title("tirefr fpc1")
    h_tirefr_fpc2=figure('Visible','off');
    scatter(domainPoints{1},fpcs(:,2),sz,'filled');
    title("tirefr fpc2")

    

    tirotb_all_fpca=FPCA(tirotb_all_cell,domainPoints,opt);
    score=getVal(tirotb_all_fpca,'xi_est');
    fpcs=getVal(tirotb_all_fpca,'phi');
    fve=getVal(tirotb_all_fpca,'FVE'); %fraction of variance explained by each pc
    pve=fve(2)*100;
    hrotb=figure('Visible','off');
    for class=1:class_num
        %rows in class
        class_rws=sum(class_count(1:class))+1:sum(class_count(1:class+1));
        %Getting first and second principle component for each
        %Coloring different classes different colors
        scatter(score(class_rws,1),score(class_rws,2),sz,colors(class,:),'filled');
        hold on
    end
    hold off
    title(append('TI by rotation, Ref ctr.',newline,string(pve), '% variance explained'));
    legend(classNames,'Interpreter','none'); %add legend for each color
    % print(hrotb,'-djpeg','tirotb_all_angiosperms_pca.jpg')
    h_tirotb_fpc1=figure('Visible','off');
    scatter(domainPoints{1},fpcs(:,1),sz,'filled');
    title("tirotb fpc1")
    h_tirotb_fpc2=figure('Visible','off');
    scatter(domainPoints{1},fpcs(:,2),sz,'filled');
    title("tirotb fpc2")

    


    tirefb_all_fpca=FPCA(tirefb_all_cell,domainPoints,opt);
    score=getVal(tirefb_all_fpca,'xi_est');
    fpcs=getVal(tirefb_all_fpca,'phi');
    fve=getVal(tirotr_all_fpca,'FVE'); %fraction of variance explained by each pc
    pve=fve(2)*100;
    hrefb=figure('Visible','off');
    for class=1:class_num
        %rows in class
        class_rws=sum(class_count(1:class))+1:sum(class_count(1:class+1));
        %Getting first and second principle component for each
        %Coloring different classes different colors
        scatter(score(class_rws,1),score(class_rws,2),sz,colors(class,:),'filled');
        hold on
    end
    hold off
    title(append('TI by reflection, Ref ctr.',newline,string(pve), '% variance explained'));
    legend(classNames,'Interpreter','none'); %add legend for each color
    % print(hrefb,'-djpeg','tirefb_all_angiosperms_pca.jpg')
    h_tirefb_fpc1=figure('Visible','off');
    scatter(domainPoints{1},fpcs(:,1),sz,'filled');
    title("tirefb fpc1")
    h_tirefb_fpc2=figure('Visible','off');
    scatter(domainPoints{1},fpcs(:,2),sz,'filled');
    title("tirefb fpc2")
    


    %Export to appropriate file
    filePath="all_pcaPlots";

    if ~exist(fullfile(pwd,filePath),"dir")
        mkdir(filePath)
    end

    save(fullfile(pwd,filePath,"data.mat"));
    
    %Data plotted against first two principal components
    exportgraphics(hrotr,fullfile(filePath,"tirotr_pca_all.jpg"));
    exportgraphics(hrefr,fullfile(filePath,"tirefr_pca_all.jpg"));
    exportgraphics(hrotb,fullfile(filePath,"tirotb_pca_all.jpg"));
    exportgraphics(hrefb,fullfile(filePath,"tirefb_pca_all.jpg"));

    savefig(hrotr,fullfile(filePath,"tirotr_pca_all.fig"));
    savefig(hrefr,fullfile(filePath,"tirefr_pca_all.fig"));
    savefig(hrotb,fullfile(filePath,"tirotb_pca_all.fig"));
    savefig(hrefb,fullfile(filePath,"tirefb_pca_all.fig"));

    %Saving the first two principal components
    exportgraphics(h_tirotr_fpc1,fullfile(filePath,"tirotr_fpc1.jpg"));
    exportgraphics(h_tirotr_fpc2,fullfile(filePath,"tirotr_fpc2.jpg"));

    exportgraphics(h_tirefr_fpc1,fullfile(filePath,"tirefr_fpc1.jpg"));
    exportgraphics(h_tirefr_fpc2,fullfile(filePath,"tirefr_fpc2.jpg"));

    exportgraphics(h_tirotb_fpc1,fullfile(filePath,"tirotb_fpc1.jpg"));
    exportgraphics(h_tirotb_fpc2,fullfile(filePath,"tirotb_fpc2.jpg"));

    exportgraphics(h_tirefb_fpc1,fullfile(filePath,"tirefb_fpc1.jpg"));
    exportgraphics(h_tirefb_fpc2,fullfile(filePath,"tirefb_fpc2.jpg"));


end












