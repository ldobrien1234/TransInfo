function ti_pca_all()


    %Loop through pcaData file
    pcaDataFolders=genpath("pcaData");
    pcaFoldersList=strsplit(pcaDataFolders,pathsep);

    %Iterating through folder paths
    for i=1:length(pcaFoldersList)
        pcaFolderPath=pcaFoldersList{i};

        if isempty(pcaFolderPath)
            continue;
        end

        pca_files=dir(fullfile(pcaFolderPath, '*.mat'));

        if isempty(pca_files)
            continue;
        end

        
        %initializing
        class_count=0; %initializing "vector" (end+1) indexing still works
        tirotr_angiosperms=[];
        tirefr_angiosperms=[];
        tirotb_angiosperms=[];
        tirefb_angiosperms=[];
        
        %Number of images in previous classes and current (in loop)
        fileNames={}; %To save filenames for legend
        for file=pca_files
            load(fullfile(pcaFolderPath,file.name),'tirotr_all','tirefr_all','tirotb_all','tirefb_all');

            %Rows are TI data for a given flower, columns are theta mesh points
            [rws,~]=size(tirotr_all); %Size of TI data function
            fileNames{end+1}=erase(file.name,'.mat');
            class_count(end+1)=rws; %Number of flowers in 1st class (early angiosperms)
            tirotr_angiosperms(end+1:end+rws,:)=tirotr_all;
            tirefr_angiosperms(end+1:end+rws,:)=tirefr_all;
            tirotb_angiosperms(end+1:end+rws,:)=tirotb_all;
            tirefb_angiosperms(end+1:end+rws,:)=tirefb_all;
            %summing previous classes with current for next loop

        end %end loop through files in a folder

    end %end loop through folders in pcaData
    
    
    class_num=length(class_count)-1; %includes empty class at beginning
    sz=50; %size of points in scatter plots
    colors=lines(class_num); %colors for each class (not including empty class)
    
    [~,score,~,~,explained] = pca(tirotr_angiosperms);
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
    title({
        ['TI by rotation, Rot ctr' ]
        [num2str(explained(1)+explained(2)) ' variance explained' ]
        });
    legend(fileNames); %add legend for each color
    % print(hrotr,'-djpeg','tirotr_all_angiosperms_pca.jpg')
    


    %Should add a legend to the plots!!
    [~,score,~,~,explained] = pca(tirefr_angiosperms);
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
    title({
        ['TI by reflection, Rot ctr' ]
        [num2str(explained(1)+explained(2)) ' variance explained' ]
        });
    legend(fileNames); %add legend for each color
    % print(hrefr,'-djpeg','tirefr_all_angiosperms_pca.jpg')
    
    [~,score,~,~,explained] = pca(tirotb_angiosperms);
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
    title({
        ['TI by rotation, Ref ctr' ]
        [num2str(explained(1)+explained(2)) ' variance explained' ]
        });
    legend(fileNames); %add legend for each color
    % print(hrotb,'-djpeg','tirotb_all_angiosperms_pca.jpg')
    
    [~,score,~,~,explained] = pca(tirefb_angiosperms);
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
    title({
        ['TI by reflection, Ref ctr'] 
        [num2str(explained(1)+explained(2)) ' variance explained' ]
        });
    legend(fileNames); %add legend for each color
    % print(hrefb,'-djpeg','tirefb_all_angiosperms_pca.jpg')

    %Export to appropriate file
    filePath="all_pcaPlots";

    if ~exist(filePath,"dir")
        mkdir(filePath)
    end


    exportgraphics(hrotr,fullfile(filePath,"tirotr_pca_all.jpg"));
    exportgraphics(hrefr,fullfile(filePath,"tirefr_pca_all.jpg"));
    exportgraphics(hrotb,fullfile(filePath,"tirotb_pca_all.jpg"));
    exportgraphics(hrefb,fullfile(filePath,"tirefb_pca_all.jpg"));



end












