function class_pca()

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
        
        [~,score,~,~,explained] = pca(tirotr_all);
        hrotr=figure('Visible', 'off'); scatter(score(:,1),score(:,2),sz,'filled')
        title({
            ['TI by rotation, Rot ctr' ] 
            [num2str(explained(1)+explained(2)) ' variance explained' ] 
            });
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirotr_all_pca.jpg')
        hrotr.UserData=filenames; %Store filenames in scatter UserData
        %Create custom data tip callback
        dcm=datacursormode(hrotr);
        dcm.Enable='on';
        dcm.UpdateFcn=@(src,event) showFilename(event,hrotr);
        
        
         [~,score,~,~,explained] = pca(tirefr_all);
        hrefr=figure('Visible', 'off'); scatter(score(:,1),score(:,2),sz,'filled')
        title({
            ['TI by reflection, Rot ctr' ] 
            [num2str(explained(1)+explained(2)) ' variance explained' ]
            });
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirefr_all_pca.jpg')
        hrefr.UserData=filenames; %Store filenames in scatter UserData
        %Create custom data tip callback
        dcm=datacursormode(hrefr);
        dcm.Enable='on';
        dcm.UpdateFcn=@(src,event) showFilename(event,hrefr);
        
        
         [~,score,~,~,explained] = pca(tirotb_all);
        hrotb=figure('Visible', 'off'); scatter(score(:,1),score(:,2),sz,'filled')
        title({
            ['TI by rotation, Ref ctr' ] 
            [num2str(explained(1)+explained(2)) ' variance explained' ]
            });
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirotb_all_pca.jpg')
        hrotb.UserData=filenames; %Store filenames in scatter UserData
        %Create custom data tip callback
        dcm=datacursormode(hrotb);
        dcm.Enable='on';
        dcm.UpdateFcn=@(src,event) showFilename(event,hrotb);
        
        
         [~,score,~,~,explained] = pca(tirefb_all);
        hrefb=figure('Visible', 'off'); scatter(score(:,1),score(:,2),sz,'filled')
        title({
            ['TI by reflection, Ref ctr' ] 
            [num2str(explained(1)+explained(2)) ' variance explained' ]
            });
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirefb_all_pca.jpg')
        hrefb.UserData=filenames; %Store filenames in scatter UserData
        %Create custom data tip callback
        dcm=datacursormode(hrefb);
        dcm.Enable='on';
        dcm.UpdateFcn=@(src,event) showFilename(event,hrefb);
        
        % tirotr_early_angiosperms=tirotr_all;
        % tirefr_early_angiosperms=tirefr_all;
        % tirotb_early_angiosperms=tirotb_all;
        % tirefb_early_angiosperms=tirefb_all;

        outputFolder="pcaData";
        plotFolder="pcaPlots";
        
        
        

        output=erase(dataFolderPath,"data"+filesep);
        [~,name,~]=fileparts(output); %Name of last folder, which becomes filename for pcaData files
        
        plotOutputPath=fullfile(plotFolder,output);
        dataOutputPath=fullfile(outputFolder,output);
        dataOutputFolder=erase(dataOutputPath,filesep+string(name)); %for creating directory


        if ~exist(dataOutputFolder, "dir")
            mkdir(dataOutputFolder);
        end

        if ~exist(plotOutputPath,"dir")
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
  


    end %end loop through folders in data folder

end

