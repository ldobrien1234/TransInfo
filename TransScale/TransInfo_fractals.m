function TransInfo_fractals(topFolder)
%% Tranformation Information - asymmetry analysis of flower images from angiosperm poster
%The sections of this matlab script go through the steps for applying TI to
%identify the approximate symmetries of the provided flower in the image. 
%It also compares to other symmetry measures based on identification of the
%center of the flower and the edges of the petals
    %topFolder is relative to where this file is stored. It should be a
    %path to a folder containing subfolders, which eventually contain
    %.png images.

    %Adding path for auxiliary functions
    parent=fileparts(pwd);
    addpath(fullfile(parent,'AuxFunctions'));


    %% Read image

    % Recursively get all folders inside topFolder (relative to pwd)
    allFolders = genpath(topFolder);  % Returns a colon-separated list of paths
    
    % Split into individual folder paths
    folderList = strsplit(allFolders, pathsep);
    
    % Loop through each folder
    for i = 1:length(folderList)
        folderPath = folderList{i};
        
        if isempty(folderPath)
            continue; %if folder is empty, pass to next iteration of loop
        end
    
        % Get png and jpg files in this folder only (non-recursive)
        pngFiles = dir(fullfile(folderPath, '*.png'));
        jpgFiles = dir(fullfile(folderPath, '*.jpg'));
        files=[pngFiles;jpgFiles]; %combine into one array
    
        % Skip folders with neither jpg's nor png's
        if isempty(files)
            continue;
        end
        
        %Keeping track of indices in loop
        file_count=0;
        
        %Loop through files in one folder
        for file=files'
        
            file_count=file_count+1;
            
            %Keeping track of path of image for storage
            imgPath=fullfile(file.folder,file.name);

            %Read and store image
            tmp=imread(imgPath);
            imdata0=tmp;
            
            
            %emphasize greenness of pixel, and add 1 to keep away from 0;
            %imdata=2*tmp(:,:,1)-tmp(:,:,2)-tmp(:,:,3)+1;
            imdata=im2gray(tmp)+1; %Convert to grayscale, take negative
            
            Dthresh=15; 
            Pmax=0;
            [M,N]=size(imdata);
            
            % % center of image 
            % pyc=(M+1)/2; 
            % pxc=(N+1)/2; 
            % 
            % %index vectors
            % %(1,1) is at the top left of the image
            % indx=1:N;
            % indy=1:M;
            % 
            % % coordinates
            % xcoord=indx-pxc; %positive right, same as indices
            % ycoord=pyc-indy; %positive up, backwards from indices
            
            
            %figure;imshow(imdata)
            
            %%%%%%%%%%%%%%%
            %% Find Center
            %%%%%%%%%%%%%%%
            
            %% reflect to find center axis of fern
            %Note: not computing TI for reflections in [0,pi] but only
            %reflecting across horizontal and vertical axes

            %y and y shifts to consider
            %Use x,y frame that starts at center of image with up and right postive
            %Only explores near the center of the image

            % ylo=ceil(min(ycoord)/4);
            % yhi=floor(max(ycoord)/4);
            % xlo=ceil(min(xcoord)/4);
            % xhi=floor(max(xcoord)/4);
            % NXmax=21; %51; %Setting grid size in x-direction
            % %NXmax=xhi-xlo+1;
            % NYmax=21; %51; %Setting grid size in y-direction
            % %NYmax=yhi-ylo+1;
            % yshifts=linspace(ylo,yhi,NYmax); 
            % xshifts=linspace(xlo,xhi,NXmax);
            % 
            % 
            % TIx=zeros(1,NXmax);
            % TIy=zeros(NYmax,1);
            % 
            % Ary=[1,0,0;0,-1,0;0,0,1]; %Affine transformation reflecting across y-axis
            % Arx=[-1,0,0;0,1,0;0,0,1]; %Affine transformation reflecting across x-axis
            % 
            % %shift in x and reflect about y
            % for llx=1:NXmax
            %         b=[xshifts(llx),0];
            %         Atb=[1,0,0;0,1,0;b(1),-b(2),1];
            %         AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
            %         Atx=AtbI*Arx*Atb;
            % 
            % 
            %         TIx(llx)=transinfo(imdata,Atx,Dthresh,Pmax);
            % 
            % end
            % 
            % for kky=1:NYmax
            %         b=[0,yshifts(kky)];
            % 
            %         Atb=[1,0,0;0,1,0;b(1),-b(2),1];
            %         AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
            %         Aty=AtbI*Ary*Atb;
            % 
            %         TIy(kky)=transinfo(imdata,Aty,Dthresh,Pmax);
            % end
            % 
            % 
            % %center will be at min of TI
            % [~,xlocs,~,proms]=findpeaks(-TIx);
            % %figure;findpeaks(-TIx,xshifts,'Annotate','extents')
            % %title('x center')
            % [~,nxc]=max(proms);
            % bxc=xshifts(xlocs(nxc)); %xshift in pixels
            % 
            % [~,ylocs,~,proms]=findpeaks(-TIy);
            % %figure;findpeaks(-TIy,yshifts,'Annotate','extents')
            % %title('y center')
            % [~,nyc]=max(proms);
            % byc=yshifts(ylocs(nyc)); %yshift in pixels

            %% Selected center
            bxc=round(N/2);
            byc=round(0.94*M);
          
            
            %% Compute TI for translation in y-direction and scaling
            % about center found by reflection

            %parameter ranges for translation and scaling
            minscale=0.5;
            scaleRange=linspace(minscale,1,30);
            transRange=linspace(0,0.8*M,30); %avoiding strange behavior near M

            TI=zeros(length(scaleRange),length(transRange));

            %Transform to center and reverse
            b=[bxc,byc];
            Atb=[1,0,0;0,1,0;-b(1),-b(2),1]; %translate to origin (center to top-left)
            AtbI=[1,0,0;0,1,0;b(1),b(2),1]; %translate from origin (top-left to center)

            for i=1:length(scaleRange)
                    s=scaleRange(i);

                for j=1:length(transRange)
                    t=transRange(j);

                    As=[s,0,0;0,s,0;0,0,1]; %scaling transformation
                    Att=[1,0,0;0,1,0;0,-t,1]; %translating in positive y-direction

                    %Total transformation (scaling and translating commute)
                    Atot=Atb*As*Att*AtbI; %Matlab applies left-to-right

                    TI(i,j)=transinfo(imdata,Atot,0,Pmax);
                end
            end
            

            
            
            h=figure('Visible', 'off'); %creating figure without displaying it
            surf(transRange,scaleRange,TI);
            xlabel('Translating');
            ylabel('Scaling');
            
            %%Writing to folder

            dataFolder='data';
            plotFolder='plots';

            if ~isfolder(dataFolder)
                mkdir(dataFolder)
            end

            if ~isfolder(plotFolder)
                mkdir(plotFolder)
            end
            
            %pwd for "path to working directory"
            
            %Gets path removing topFolder/ (e.g. angiosperms/)
            relativeFolder=erase(folderPath, topFolder);

            %Create output file path
            dataOutput=fullfile(dataFolder,relativeFolder);
            plotOutput=fullfile(plotFolder,relativeFolder);

            %Constructing output folders
            if ~exist(dataOutput,'dir')
                mkdir(dataOutput)
            end

            if ~exist(plotOutput, 'dir')
                mkdir(plotOutput)
            end
            
            
            %Naming files and giving them paths
            [~,name,~]=fileparts(file.name);
            
            plotTitle=[name,'_plots.jpg'];
            dataTitle=[name,'_data.mat'];

            plotFilePath=fullfile(plotOutput,plotTitle);
            dataFilePath=fullfile(dataOutput,dataTitle);

            %Writing files to appropriate path
            exportgraphics(h,plotFilePath); %image file
            save(dataFilePath); %saving matlab workspace to file path

            close(h);
            
            
        end %end of png files loop
    end %end of folders loop
    
end

        
        % ax=subplot(1,3,2);
        % hold on;
        % image(imdata),colormap(ax,flipud(gray(256))), 
        % axis image, axis off; 
        % set(ax,'YDir','reverse','DataAspectRatio',[1 1 1]);
        % xlim(xClim); ylim(yClim);
        % 
        %  
        % 
        % line(pxc*[1,1], yClim,'Color','k','Linewidth',lw,'LineStyle',':')
        % line(xClim,pyc*[1 1],'Color','k','Linewidth',lw,'LineStyle',':')
        % %rotation center
        % line(xClim,(pyc-ryc)*[1 1],'Color',lrcolor,'Linewidth',lw/2,'LineStyle','--')
        % plot((pxc+rxc),(pyc-ryc),'Color',lrcolor,'Marker','*','MarkerSize',ms)
        % plot((pxc+rxc),(pyc-ryc),'Color',lrcolor,'Marker','o','MarkerSize',ms)
        % %translation center
        % line(xClim,(pyc-byc)*[1 1],'Color',ltcolor,'Linewidth',lw/2,'LineStyle','--')
        % plot((pxc+bxc),(pyc-byc),'Color',ltcolor,'Marker','x','MarkerSize',ms)
        % plot((pxc+bxc),(pyc-byc),'Color',ltcolor,'Marker','o','MarkerSize',ms)
        % %area center
        % plot((pxc+axc),(pyc-ayc),'Color','k','Marker','^','MarkerSize',ms/2)
        % plot((pxc+axc),(pyc-ayc),'Color','k','Marker','o','MarkerSize',ms)
        % %manual center
        % plot((pxc+xFc),(pyc-yFc),'Color','k','Marker','.','MarkerSize',ms)
        % plot((pxc+xFc),(pyc-yFc),'Color','k','Marker','o','MarkerSize',ms)