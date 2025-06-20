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
    
        for file=data_files_all'
           load(fullfile(dataFolderPath,file.name), "TIrotr","TIrefr","TIrotb","TIrefb");
           tirotr_all(file_count_all,:)=TIrotr;
           tirefr_all(file_count_all,:)=TIrefr;
           tirotb_all(file_count_all,:)=TIrotb;
           tirefb_all(file_count_all,:)=TIrefb;
           file_count_all=file_count_all+1;
        end
        
        %c = linspace(1,10,file_count-1);
        c = jet(file_count_all-1);
        sz=50;
        
        [~,score,~,~,explained] = pca(tirotr_all);
        hrotr=figure('Visible', 'off'); scatter(score(:,1),score(:,2),sz,c,'filled')
        title({
            ['TI by rotation, Rot ctr' ] 
            [num2str(explained(1)+explained(2)) ' variance explained' ] 
            });
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirotr_all_pca.jpg')
        
         [~,score,~,~,explained] = pca(tirefr_all);
        hrefr=figure('Visible', 'off'); scatter(score(:,1),score(:,2),sz,c,'filled')
        title({
            ['TI by reflection, Rot ctr' ] 
            [num2str(explained(1)+explained(2)) ' variance explained' ]
            });
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirefr_all_pca.jpg')
        
        
         [~,score,~,~,explained] = pca(tirotb_all);
        hrotb=figure('Visible', 'off'); scatter(score(:,1),score(:,2),sz,c,'filled')
        title({
            ['TI by rotation, Ref ctr' ] 
            [num2str(explained(1)+explained(2)) ' variance explained' ]
            });
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirotb_all_pca.jpg')
        
        
         [~,score,~,~,explained] = pca(tirefb_all);
        hrefb=figure('Visible', 'off'); scatter(score(:,1),score(:,2),sz,c,'filled')
        title({
            ['TI by reflection, Ref ctr' ] 
            [num2str(explained(1)+explained(2)) ' variance explained' ]
            });
        %ftitle=[file.name,'tirotr_all_pca.jpg'];
        % print(h,'-djpeg','tirefb_all_pca.jpg')
        
        % tirotr_early_angiosperms=tirotr_all;
        % tirefr_early_angiosperms=tirefr_all;
        % tirotb_early_angiosperms=tirotb_all;
        % tirefb_early_angiosperms=tirefb_all;

        outputFolder="pcaData";
        plotFolder="pcaPlots";
        
        
        

        output=erase(dataFolderPath,["data",filesep]);
        [~,name,~]=fileparts(output); %Name of last folder, which becomes filename for pcaData files
        
        plotOutputPath=fullfile(plotFolder,output);
        dataOutputPath=fullfile(outputFolder,output);
        dataOutputFolder=erase(dataOutputPath,[filesep,name]); %for creating directory


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
     



    end %end loop through folders in data folder

end

