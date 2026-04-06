function class_pca()


    %Path for functional PCA package
    pacePath="/Users/liamobrien/Documents/PACE_matlab-master/release2.17/PACE";
    addpath(genpath(pacePath));
    %Option to get 2 principal components each time
    opt = setOptions();
    opt.selection_k = 'FVE';
    opt.maxk = 2;
    opt.FVE_threshold = 0.99;
    


    %Getting paths from data folder
    allDataFolders=genpath("data");
    dataFolderList = strsplit(allDataFolders, pathsep);

    %Iterating through folder paths)
    for i = 1:length(dataFolderList)
        dataFolderPath = dataFolderList{i};
        
        if isempty(dataFolderPath)
            continue; %continue passes to next iteration of loop
        end
    
        % Get data files in this folder only (non-recursive)
        data_files_all = dir(fullfile(dataFolderPath, '*data.mat'));

        if isempty(data_files_all)
            continue;
        end
    
        %files=dir('15.png');
        file_count_all=1;

        tirotr_all=[];
        tirefr_all=[];
        tirotb_all=[];
        tirefb_all=[];
        
        filenames={};
        for file=data_files_all'
           filenames{end+1}=erase(file.name,'_data.mat');
           load(fullfile(dataFolderPath,file.name), "TIrotr","TIrefr","TIrotb","TIrefb");
           tirotr_all(file_count_all,:)=TIrotr;
           tirefr_all(file_count_all,:)=TIrefr;
           tirotb_all(file_count_all,:)=TIrotb;
           tirefb_all(file_count_all,:)=TIrefb;
           file_count_all=file_count_all+1;
        end
        
        %c = linspace(1,10,file_count-1);
        %c = jet(file_count_all-1);
        sz=50;

        %Converting matrices to cell arrays for functional pca
        tirotr_all_cell=num2cell(tirotr_all,2)';
        tirefr_all_cell=num2cell(tirefr_all,2)';
        tirotb_all_cell=num2cell(tirotb_all,2)';
        tirefb_all_cell=num2cell(tirefb_all,2)';

        domainPoints=linspace(0,2*pi,length(tirotr_all(1,:)));
        %convert to cell that matches the number of TI curves
        domainPoints=repmat({domainPoints}, numel(tirotr_all_cell), 1)';
        
        tirotr_all_fpca=FPCA(tirotr_all_cell,domainPoints,opt);
        score=getVal(tirotr_all_fpca,'xi_est');
        fpcs=getVal(tirotr_all_fpca,'phi');
        fve=getVal(tirotr_all_fpca,'FVE'); %fraction of variance explained by each pc
        pve=fve(2)*100;
        hrotr=figure('Visible','off'); 
        scatter(score(:,1),score(:,2),sz,'filled')
        title(append('TI by rotation, Rot ctr.',newline,string(pve),'% variance explained')); 
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirotr_all_pca.jpg')
        hrotr.UserData=filenames; %Store filenames in scatter UserData
        %Create custom data tip callback
        dcm=datacursormode(hrotr);
        dcm.Enable='on';
        dcm.UpdateFcn=@(src,event) showFilename(event,hrotr);
        hrotr_fpc1=figure('Visible','Off');
        title("tirotr fpc1")
        scatter(domainPoints{1},fpcs(:,1),sz,'filled');
        hrotr_fpc2=figure('Visible','Off');
        scatter(domainPoints{1},fpcs(:,2),sz,'filled');
        title("tirotr fpc2")
        
        
        tirefr_all_fpca=FPCA(tirefr_all_cell,domainPoints,opt);
        score=getVal(tirefr_all_fpca,'xi_est');
        fpcs=getVal(tirefr_all_fpca,'phi'); 
        fve=getVal(tirefr_all_fpca,'FVE'); %fraction of variance explained by each pc
        pve=fve(2)*100;
        hrefr=figure('Visible', 'off'); 
        scatter(score(:,1),score(:,2),sz,'filled')
        title(append('TI by reflection, Rot ctr.',newline,string(pve),'% variance explained')); 
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirefr_all_pca.jpg')
        hrefr.UserData=filenames; %Store filenames in scatter UserData
        %Create custom data tip callback
        dcm=datacursormode(hrefr);
        dcm.Enable='on';
        dcm.UpdateFcn=@(src,event) showFilename(event,hrefr);
        hrefr_fpc1=figure('Visible','Off');
        title("tirefr fpc1")
        scatter(domainPoints{1},fpcs(:,1),sz,'filled');
        hrefr_fpc2=figure('Visible','Off');
        scatter(domainPoints{1},fpcs(:,2),sz,'filled');
        title("tirefr fpc2")
        
        
        tirotb_all_fpca=FPCA(tirotb_all_cell,domainPoints,opt);
        score=getVal(tirotb_all_fpca,'xi_est');
        fpcs=getVal(tirotb_all_fpca,'phi');
        fve=getVal(tirotb_all_fpca,'FVE'); %fraction of variance explained by each pc
        pve=fve(2)*100;
        hrotb=figure('Visible', 'off');
        scatter(score(:,1),score(:,2),sz,'filled')
        title(append('TI by rotation, Ref ctr.',newline,string(pve),'% variance explained'));
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirotb_all_pca.jpg')
        hrotb.UserData=filenames; %Store filenames in scatter UserData
        %Create custom data tip callback
        dcm=datacursormode(hrotb);
        dcm.Enable='on';
        dcm.UpdateFcn=@(src,event) showFilename(event,hrotb);
        hrotb_fpc1=figure('Visible','Off');
        title("tirotb fpc1")
        scatter(domainPoints{1},fpcs(:,1),sz,'filled');
        hrotb_fpc2=figure('Visible','Off');
        scatter(domainPoints{1},fpcs(:,2),sz,'filled');
        title("tirotb fpc2")
        
        
        tirefb_all_fpca=FPCA(tirefb_all_cell,domainPoints,opt);
        score=getVal(tirefb_all_fpca,'xi_est');
        fpcs=getVal(tirefb_all_fpca,'phi');
        fve=getVal(tirefb_all_fpca,'FVE'); %fraction of variance explained by each pc
        pve=fve(2)*100;        
        hrefb=figure('Visible', 'off');
        scatter(score(:,1),score(:,2),sz,'filled')
        title(append('TI by reflection, Ref ctr.',newline,string(pve),'% variance explained'));
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirefb_all_pca.jpg')
        hrefb.UserData=filenames; %Store filenames in scatter UserData
        %Create custom data tip callback
        dcm=datacursormode(hrefb);
        dcm.Enable='on';
        dcm.UpdateFcn=@(src,event) showFilename(event,hrefb);
        hrefb_fpc1=figure('Visible','Off');
        title("tirefb fpc1")
        scatter(domainPoints{1},fpcs(:,1),sz,'filled');
        hrefb_fpc2=figure('Visible','Off');
        scatter(domainPoints{1},fpcs(:,2),sz,'filled');
        title("tirefb fpc2")
        
        % tirotr_early_angiosperms=tirotr_all;
        % tirefr_early_angiosperms=tirefr_all;
        % tirotb_early_angiosperms=tirotb_all;
        % tirefb_early_angiosperms=tirefb_all;

        outputFolder="pcaData";
        plotFolder="pcaPlots";
        
        
        

        output=erase(dataFolderPath,"data"+filesep);
        [~,name,~]=fileparts(output); %Name of last folder, which becomes filename for pcaData files
        
        plotOutputPath=fullfile(pwd,plotFolder,output);
        dataOutputPath=fullfile(pwd,outputFolder,output);
        dataOutputFolder=erase(dataOutputPath,filesep+string(name)); %for creating directory


        if ~exist(fullfile(pwd,dataOutputFolder), "dir")
            mkdir(dataOutputFolder);
        end

        if ~exist(fullfile(pwd,plotOutputPath),"dir")
            mkdir(plotOutputPath)
        end


        save(dataOutputPath);
        exportgraphics(hrotr,fullfile(plotOutputPath,"tirotr_pca.jpg"));
        exportgraphics(hrefr,fullfile(plotOutputPath,"tirefr_pca.jpg"));
        exportgraphics(hrotb,fullfile(plotOutputPath,"tirotb_pca.jpg"));
        exportgraphics(hrefb,fullfile(plotOutputPath,"tirefb_pca.jpg"));
        
        %Saving interactive matlab figures
        savefig(hrotr, fullfile(plotOutputPath,"tirotr_pca.fig"));
        savefig(hrefr, fullfile(plotOutputPath,"tirefr_pca.fig"));
        savefig(hrotb, fullfile(plotOutputPath,"tirotb_pca.fig"));
        savefig(hrefb, fullfile(plotOutputPath,"tirefb_pca.fig"));

            %Saving the first two principal components
        exportgraphics(hrotr_fpc1,fullfile(plotOutputPath,"tirotr_fpc1.jpg"));
        exportgraphics(hrotr_fpc2,fullfile(plotOutputPath,"tirotr_fpc2.jpg"));

        exportgraphics(hrefr_fpc1,fullfile(plotOutputPath,"tirefr_fpc1.jpg"));
        exportgraphics(hrefr_fpc2,fullfile(plotOutputPath,"tirefr_fpc2.jpg"));

        exportgraphics(hrotb_fpc1,fullfile(plotOutputPath,"tirotb_fpc1.jpg"));
        exportgraphics(hrotb_fpc2,fullfile(plotOutputPath,"tirotb_fpc2.jpg"));

        exportgraphics(hrefb_fpc1,fullfile(plotOutputPath,"tirefb_fpc1.jpg"));
        exportgraphics(hrefb_fpc2,fullfile(plotOutputPath,"tirefb_fpc2.jpg"));
  


    end %end loop through folders in data folder

end

