function ti_pca_all()


    %Loop through pcaData file
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
    
    
    class_num=length(class_count)-1; %includes empty class at beginning
    sz=50; %size of points in scatter plots
    colors=lines(class_num); %colors for each class (not including empty class)
    
    [~,score,~,~,explained] = pca(tirotr_all_data);
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
    legend(classNames,'Interpreter','none'); %add legend for each color
    % print(hrotr,'-djpeg','tirotr_all_angiosperms_pca.jpg')
    


    %Should add a legend to the plots!!
    [~,score,~,~,explained] = pca(tirefr_all_data);
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
    legend(classNames,'Interpreter','none'); %add legend for each color
    % print(hrefr,'-djpeg','tirefr_all_angiosperms_pca.jpg')
    
    [~,score,~,~,explained] = pca(tirotb_all_data);
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
    legend(classNames,'Interpreter','none'); %add legend for each color
    % print(hrotb,'-djpeg','tirotb_all_angiosperms_pca.jpg')
    
    [~,score,~,~,explained] = pca(tirefb_all_data);
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
    legend(classNames,'Interpreter','none'); %add legend for each color
    % print(hrefb,'-djpeg','tirefb_all_angiosperms_pca.jpg')

    %Export to appropriate file
    filePath="all_pcaPlots";

    if ~exist(filePath,"dir")
        mkdir(filePath)
    end

    save(fullfile(filePath,"data.mat"));

    exportgraphics(hrotr,fullfile(filePath,"tirotr_pca_all.jpg"));
    exportgraphics(hrefr,fullfile(filePath,"tirefr_pca_all.jpg"));
    exportgraphics(hrotb,fullfile(filePath,"tirotb_pca_all.jpg"));
    exportgraphics(hrefb,fullfile(filePath,"tirefb_pca_all.jpg"));

    savefig(hrotr,fullfile(filePath,"tirotr_pca_all.fig"));
    savefig(hrefr,fullfile(filePath,"tirefr_pca_all.fig"));
    savefig(hrotb,fullfile(filePath,"tirotb_pca_all.fig"));
    savefig(hrefb,fullfile(filePath,"tirefb_pca_all.fig"));



end












