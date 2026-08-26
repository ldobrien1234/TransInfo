function TransInfo_images(topFolder)
%% Tranformation Information - asymmetry analysis of flower images from angiosperm poster
%The sections of this matlab script go through the steps for applying TI to
%identify the approximate symmetries of the provided flower in the image. 
%It also compares to other symmetry measures based on identification of the
%center of the flower and the edges of the petals.
    %topFolder is relative to where this file is stored. It should be a
    %path to a folder containing subfolders, which eventually contain
    %.png images.

    %Adding path for auxiliary functions
    parent=fileparts(pwd);
    addpath(fullfile(parent,'AuxFunctions'));


    %% Read image

    % Recursively get all folders inside topFolder
    parent=fileparts(pwd);
    fullTop=fullfile(parent,topFolder);
    allFolders = genpath(fullTop);  % Returns a colon-separated list of paths
    
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
    
        % Skip folders with no PNGs
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
            imdata=im2gray(tmp)+1; %Convert to grayscale
            
            Dthresh=15; 
            Pmax=256;
            [M,N]=size(imdata);
            
            % center of image 
            pyc=ceil(M/2); 
            pxc=ceil(N/2); 
            
            %index vectors
            %(1,1) is at the top left of the image 
            % coordinates to search through rectangle in the middle of the image
            xcoord=(pxc-floor(N/4)):(pxc+ceil(N/4)); %positive right, same as indices
            ycoord=(pyc-floor(M/4)):(pyc+ceil(M/4)); %positive up, backwards from indices
            
            
            %figure;imshow(imdata)
            
            %%%%%%%%%%%%%%%
            %% Find Center
            %%%%%%%%%%%%%%%
            
            %% Find center by rotation about x,y
            
            %y and y shifts to consider
            %Use x,y frame that starts at center of image with up and right postive
            ylo=min(ycoord);
            yhi=max(ycoord);
            xlo=min(xcoord);
            xhi=max(xcoord);
            NXmax=21; %51; %Setting grid size in x-direction
            %NXmax=xhi-xlo+1;
            NYmax=21; %51; %Setting grid size in y-direction
            %NYmax=yhi-ylo+1;
            yshifts=round(linspace(ylo,yhi,NYmax)); 
            xshifts=round(linspace(xlo,xhi,NXmax));
            
            %angles of rotation to consider
            dTh=pi/6; 
            thetavals=dTh:dTh:(2*pi-dTh); %Measure angle conterclockwise from X
            
            %downsample image (i.e. coarsen mesh)
            %npskip=16;
            npskip=1;
            
            b=zeros(1,2);

            %Stores TI for all centers in the grid
            TIcenter=zeros(NXmax,NYmax,length(thetavals));
            COVcenter=zeros(NXmax,NYmax,length(thetavals));

            %Candidates whose transformed domain overlaps less than
            %COVthresh of the original flower's domain are rejected:
            %transinfo's TI is an *average* over that overlap, so a
            %wrong-but-far center that rotates most of the flower out of
            %its own footprint can leave only a handful of leftover
            %pixels whose average divergence is spuriously small,
            %mimicking a true symmetry center. See transinfo.m's
            %coverage output for details.
            COVthresh=0.7;

            for llx=1:NXmax
                b(1)=xshifts(llx)/npskip;
                for kky=1:NYmax
                    b(2)=yshifts(kky)/npskip;

                    %shift right and up
                    Atb=[1,0,0;0,1,0;-b(1),-b(2),1];
                    AtbI=[1,0,0;0,1,0;b(1),b(2),1];

                    for mth=1:length(thetavals)
                        theta=thetavals(mth);
                        %rotate theta ccw about origin
                        Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
                        %ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];


                        Atform=Atb*Arot*AtbI; %transforms right multiplied, so leftmost is carried out first
                        [TIcenter(llx,kky,mth),COVcenter(llx,kky,mth)]=transinfo(imdata(1:npskip:end,1:npskip:end),Atform,Dthresh,Pmax);
                    end

                end

            end

            %figure;imagesc(xshifts,yshifts,TIcenter(:,:,1)')
            %figure;imagesc(xshifts,yshifts,(max(TIcenter,[],3)-min(TIcenter,[],3))')
            %figure;imagesc(xshifts,yshifts,sum(abs(diff(TIcenter,[],3)).^2,3)')

            %arr=sum(abs(diff(TIcenter,[],3)).^2,3)';

            TIcenterMasked=TIcenter;
            TIcenterMasked(COVcenter<COVthresh)=Inf;

            %For each center point in the grid, find the minimum
            for rarr=1:NYmax
                for carr=1:NXmax
                    arr(rarr,carr)=min(TIcenterMasked(rarr,carr,:));
                end
            end

            %arr=TIcenter;
            minimum=min(min(arr)); %Find the minimum of the mins
            if ~isfinite(minimum)
                %No candidate anywhere in the grid kept enough of the
                %flower's domain; fall back to the unmasked search
                %rather than erroring out.
                warning('TransInfo_images:noCoverage', ...
                    'No candidate center reached %.0f%% domain coverage; falling back to unmasked TI.', COVthresh*100);
                for rarr=1:NYmax
                    for carr=1:NXmax
                        arr(rarr,carr)=min(TIcenter(rarr,carr,:));
                    end
                end
                minimum=min(min(arr));
            end
            [xidmax,yidmax]=find(arr==minimum,1); %index of center
            %xidmax=mod(xidmaxtmp,NXmax);
            rxc=xshifts(xidmax); %center value in relation to midpoint
            ryc=yshifts(yidmax);
            %rxc=xshifts(11);
            %ryc=yshifts(13);
            
            %% translate and reflect x,y to find center of TI
            %Note: not computing TI for reflections in [0,pi] but only
            %reflecting across horizontal and vertical axes
            
            
            TIx=zeros(1,NXmax);
            TIy=zeros(NYmax,1);
            COVx=zeros(1,NXmax);
            COVy=zeros(NYmax,1);

            Ary=[1,0,0;0,-1,0;0,0,1]; %Affine transformation reflecting across x-axis (y mapsto -y)
            Arx=[-1,0,0;0,1,0;0,0,1]; %Affine transformation reflecting across y-axis (x mapsto -x)

            %shift in x and reflect across y-axis
            for llx=1:NXmax
                    b=[xshifts(llx),0];
                    Atb=[1,0,0;0,1,0;-b(1),-b(2),1];
                    AtbI=[1,0,0;0,1,0;b(1),b(2),1];
                    Atx=Atb*Arx*AtbI;


                    [TIx(llx),COVx(llx)]=transinfo(imdata,Atx,Dthresh,Pmax);

            end

            %shift in y and reflect across x-axis
            for kky=1:NYmax
                    b=[0,yshifts(kky)];

                    Atb=[1,0,0;0,1,0;-b(1),-b(2),1];
                    AtbI=[1,0,0;0,1,0;b(1),b(2),1];
                    Aty=Atb*Ary*AtbI;

                    [TIy(kky),COVy(kky)]=transinfo(imdata,Aty,Dthresh,Pmax);
            end


            %center will be at min of TI. Reject candidates below
            %COVthresh (defined above, in the rotation-center search) the
            %same way, falling back to the unmasked TIx/TIy if none
            %qualify.
            TIxMasked=TIx;
            TIxMasked(COVx<COVthresh)=Inf;
            if ~isfinite(min(TIxMasked))
                warning('TransInfo_images:noCoverage', ...
                    'No x-reflection candidate reached %.0f%% domain coverage; falling back to unmasked TI.', COVthresh*100);
                TIxMasked=TIx;
            end
            [pks,xlocs,~,~]=findpeaks(-TIxMasked);
            %Pad peaks and indices to include endpoints
            pks=horzcat(TIxMasked(1),-pks,TIxMasked(end));
            xlocs=horzcat(1,xlocs,length(-TIxMasked));
            %figure;findpeaks(-TIx,xshifts,'Annotate','extents')
            %title('x center')
            [~,nxc]=min(pks);
            bxc=xshifts(xlocs(nxc)); %xshift in pixels

            TIyMasked=TIy;
            TIyMasked(COVy<COVthresh)=Inf;
            if ~isfinite(min(TIyMasked))
                warning('TransInfo_images:noCoverage', ...
                    'No y-reflection candidate reached %.0f%% domain coverage; falling back to unmasked TI.', COVthresh*100);
                TIyMasked=TIy;
            end
            [pks,ylocs,~,~]=findpeaks(-TIyMasked);
            pks=pks';
            ylocs=ylocs';
            %Pad peaks and indices to include endpoints
            pks=horzcat(TIyMasked(1),-pks,TIyMasked(end));
            ylocs=horzcat(1,ylocs,length(-TIyMasked));
            %figure;findpeaks(-TIy,yshifts,'Annotate','extents')
            %title('y center')
            [~,nyc]=min(pks);
            byc=yshifts(ylocs(nyc)); %yshift in pixels
            
            %get min and max y locations of peaks on plateau
            %bycMin=min(yshifts(ylocs));
            %bycMax=max(yshifts(ylocs));
            %byc2=120; %identify from graph
            
            %figure; imagesc(xcoord,ycoord,imdata)
            %colormap(gray)
            %hold on
            %plot(rxc,ryc,'ro')
            %plot(bxc,byc,'c*')
            %hold off
            
            %% rotate and reflect about center found by rotation
            
            %b=[bxc,byc];
            b=[rxc,ryc];
            %b=[xFc,yFc]; %geometric centre of flower
            Atb=[1,0,0;0,1,0;-b(1),-b(2),1];
            AtbI=[1,0,0;0,1,0;b(1),b(2),1];
            
            Arx=[-1,0,0;0,1,0;0,0,1]; %reflect across y-axis
            Ary=[1,0,0;0,-1,0;0,0,1]; %reflect across x-axis
            
            Nthmax=361;
            thetavals=linspace(0,360,Nthmax)*pi/180;
            
            
            TIrotr=zeros(Nthmax,1);
            TIrefr=zeros(Nthmax,1);
            
            for llth=1:Nthmax
            
                theta=thetavals(llth);
                Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
                ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
            
                
                %rotate by theta
                tform=Atb*Arot*AtbI;
                TIrotr(llth)=transinfo(imdata,tform,Dthresh,Pmax);
                
                %reflect about axis at theta
                tform=Atb*ArotI*Ary*Arot*AtbI; %bring axis to x-axis, reflect y->-y, rotate back
                TIrefr(llth)=transinfo(imdata,tform,Dthresh,Pmax);
                
               
                
                 
                
             end
                
            % figure;
            % subplot(1,2,1),polarplot(thetavals,TIrotr'),title('rotation');
            % subplot(1,2,2),polarplot(thetavals,TIrefr'),title('reflection');
            % 
            % 
            % figure;
            % subplot(2,1,1),plot(thetavals*180/pi,TIrotr'),title('rotation');
            % subplot(2,1,2),plot(thetavals*180/pi,TIrefr'),title('reflection');
            
            %% rotate and reflect about center found by reflection
            
            b=[bxc,byc];
            %b=[rxc,ryc];
            %b=[xFc,yFc]; %geometric centre of flower
            Atb=[1,0,0;0,1,0;-b(1),-b(2),1];
            AtbI=[1,0,0;0,1,0;b(1),b(2),1];
            
            Arx=[-1,0,0;0,1,0;0,0,1];
            Ary=[1,0,0;0,-1,0;0,0,1];
            
            Nthmax=361;
            thetavals=linspace(0,360,Nthmax)*pi/180;
            
            
            TIrotb=zeros(Nthmax,1);
            TIrefb=zeros(Nthmax,1);
            
            for llth=1:Nthmax
            
                theta=thetavals(llth);
                Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
                ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
            
                
                %rotate by theta
                tform=Atb*Arot*AtbI;
                TIrotb(llth)=transinfo(imdata,tform,Dthresh,Pmax);
                
                %reflect about axis at theta
                tform=Atb*ArotI*Ary*Arot*AtbI; %bring axis to x-axis, reflect y->-y, rotate back
                TIrefb(llth)=transinfo(imdata,tform,Dthresh,Pmax);
                
               
                
                 
                
             end
                
            % figure;
            % subplot(1,2,1),polarplot(thetavals,TIrotb'),title('rotation');
            % subplot(1,2,2),polarplot(thetavals,TIrefb'),title('reflection');
            
            
            
            h=figure('Visible', 'off'); %creating figure without displaying it
            subplot(3,2,1), imshow(imdata0);
            %subplot(3,2,2), imagesc(xcoord,ycoord,imdata), colormap(gray), axis(equal), hold on, plot(rxc,ryc,'ro'), plot(bxc,byc,'c*'), hold off;
            ax=subplot(3,2,2); hold on, image(imdata), colormap(gray); axis image, axis off; 
            set(ax,'YDir','reverse','DataAspectRatio',[1 1 1]);
            plot(rxc,ryc,'ro');
            plot(bxc,byc,'c*');
            hold off
            subplot(3,2,3),plot(thetavals*180/pi,TIrotr'),title('rotation - ctr rot (red)');
            subplot(3,2,4),plot(thetavals*180/pi,TIrefr'),title('reflection - ctr rot (red)');
            subplot(3,2,5),plot(thetavals*180/pi,TIrotb'),title('rotation - ctr ref (cyan)');
            subplot(3,2,6),plot(thetavals*180/pi,TIrefb'),title('reflection - ctr ref (cyan)');
            %print('FlowerSymmetries.png','-dpng')

            %%Writing to folder

            dataFolder='data';
            plotFolder='plots';

            if exist(fullfile(pwd,dataFolder),'dir')==0
                mkdir(dataFolder)
            end

            if exist(fullfile(pwd,plotFolder),'dir')==0
                mkdir(plotFolder)
            end
            
            %pwd for "path to working directory"
            
            %Gets path removing topFolder/ (e.g. angiosperms/)
            relativeFolder=erase(folderPath, fullTop);

            %Create output file path
            dataOutput=fullfile(dataFolder,relativeFolder);
            plotOutput=fullfile(plotFolder,relativeFolder);

            %Constructing output folders
            if ~exist(fullfile(pwd,dataOutput),'dir')
                mkdir(dataOutput)
            end

            if ~exist(fullfile(pwd,plotOutput), 'dir')
                mkdir(plotOutput)
            end
            
            
            %Naming files and giving them paths
            [~,name,~]=fileparts(file.name);
            
            plotTitle=[name,'_plots.jpg'];
            dataTitle=[name,'_data.mat'];

            plotFilePath=fullfile(pwd,plotOutput,plotTitle);
            dataFilePath=fullfile(pwd,dataOutput,dataTitle);

            %Writing files to appropriate path
            exportgraphics(h,plotFilePath); %image file
            save(dataFilePath); %saving matlab workspace to file path

            close(h);
            
            
        end %end of png files loop
    end %end of folders loop
    
end

        