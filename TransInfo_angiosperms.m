function TransInfo_angiosperms(topFolder)
%% Tranformation Information - asymmetry analysis of flower images from angiosperm poster
%The sections of this matlab script go through the steps for applying TI to
%identify the approximate symmetries of the provided flower in the image. 
%It also compares to other symmetry measures based on identification of the
%center of the flower and the edges of the petals.
    %topFolder is relative to where this file is stored. It should be a
    %path to a folder containing subfolders, which eventually contain
    %.png images.




    %% Read image

    % Recursively get all folders inside topFolder
    allFolders = genpath(topFolder);  % Returns a colon-separated list of paths
    
    % Split into individual folder paths
    folderList = strsplit(allFolders, pathsep);
    
    % Loop through each folder
    for i = 1:length(folderList)
        folderPath = folderList{i};
        
        if isempty(folderPath)
            continue; %continue passes to next iteration of loop
        end
    
        % Get png and jpg files in this folder only (non-recursive)
        pngFiles = dir(fullfile(folderPath, '*.png'));
        jpgFiles = dir(fullfile(folderPath, '*.jpg'));
        files=[pngFiles;jpgFiles]; %combine into one array
    
        % Skip folders with no PNGs
        if isempty(files)
            continue;
        end

        file_count=0;
        
        for file=files'
        
            file_count=file_count+1;
            
            imgPath=fullfile(file.folder,file.name);
            tmp=imread(imgPath);
            imdata0=tmp;
            
            
            %emphasize greenness of pixel, and add 1 to keep away from 0;
            %imdata=2*tmp(:,:,1)-tmp(:,:,2)-tmp(:,:,3)+1;
            imdata=rgb2gray(tmp)+1;
            
            Dthresh=1;
            Pmax=256;
            [M,N]=size(imdata);
            
            % center of image 
            pyc=(M+1)/2; 
            pxc=(N+1)/2; 
            
            %index vectors
            %(1,1) is at the top left of the image
            indx=1:N;
            indy=1:M;
            
            % coordinates
            xcoord=indx-pxc; %positive right, same as indices
            ycoord=pyc-indy; %positive up, backwards from indices
            
            
            %figure;imshow(imdata)
            
            %%%%%%%%%%%%%%%
            %% Find Center
            %%%%%%%%%%%%%%%
            
            %% Find center by rotation about x,y
            
            %y and y shifts to consider
            %Use x,y frame that starts at center of image with up and right postive
            %Only explores near the center of the image
            ylo=ceil(min(ycoord)/4);
            yhi=floor(max(ycoord)/4);
            xlo=ceil(min(xcoord)/4);
            xhi=floor(max(xcoord)/4);
            NXmax=21; %51;
            %NXmax=xhi-xlo+1;
            NYmax=21; %51;
            %NYmax=yhi-ylo+1;
            yshifts=linspace(ylo,yhi,NYmax); 
            xshifts=linspace(xlo,xhi,NXmax);
            
            %angles of rotation to consider
            dTh=2*pi/12; 
            thetavals=dTh:dTh:(2*pi-dTh); %Measure angle conterclockwise from X
            
            %downsample image (i.e. coarsen mesh)
            %npskip=16;
            npskip=1;
            
            b=zeros(1,2);
            TIcenter=zeros(NXmax,NYmax,length(thetavals));
            
            
            for llx=1:NXmax
                b(1)=xshifts(llx)/npskip;
                for kky=1:NYmax
                    b(2)=yshifts(kky)/npskip;
                    
                    %shift right and up
                    Atb=[1,0,0;0,1,0;b(1),-b(2),1];
                    AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
                    
                    for mth=1:length(thetavals)
                        theta=thetavals(mth);
                        %rotate theta ccw about origin
                        Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
                        %ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
            
                        
                        Atform=AtbI*Arot*Atb; %transforms right multiplied, so leftmost is carried out first
                        TIcenter(llx,kky,mth)=transinfo(imdata(1:npskip:end,1:npskip:end),Atform,Dthresh,Pmax);
                    end
                    
                end
                
            end
            
            %figure;imagesc(xshifts,yshifts,TIcenter(:,:,1)')
            %figure;imagesc(xshifts,yshifts,(max(TIcenter,[],3)-min(TIcenter,[],3))')
            %figure;imagesc(xshifts,yshifts,sum(abs(diff(TIcenter,[],3)).^2,3)')
            
            %arr=sum(abs(diff(TIcenter,[],3)).^2,3)';
            for rarr=1:NYmax
                for carr=1:NXmax
                    arr(rarr,carr)=min(TIcenter(rarr,carr,:));
                end
            end
            
            %arr=TIcenter;
            minimum=min(min(arr));
            [xidmax,yidmax]=find(arr==minimum);
            %xidmax=mod(xidmaxtmp,NXmax);
            rxc=xshifts(xidmax);
            ryc=yshifts(yidmax);
            %rxc=xshifts(11);
            %ryc=yshifts(13);
            
            %% translate and reflect x,y to find center by TI
            
            
            TIx=zeros(1,NXmax);
            TIy=zeros(NYmax,1);
            
            Ary=[1,0,0;0,-1,0;0,0,1];
            Arx=[-1,0,0;0,1,0;0,0,1];
            
            %shift in x and reflect about y
            for llx=1:NXmax
                    b=[xshifts(llx),0];
                    Atb=[1,0,0;0,1,0;b(1),-b(2),1];
                    AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
                    Atx=AtbI*Arx*Atb;
                    
                    
                    TIx(llx)=transinfo(imdata,Atx,Dthresh,Pmax);
             
            end
            
            for kky=1:NYmax
                    b=[0,yshifts(kky)];
                    
                    Atb=[1,0,0;0,1,0;b(1),-b(2),1];
                    AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
                    Aty=AtbI*Ary*Atb;
                    
                    TIy(kky)=transinfo(imdata,Aty,Dthresh,Pmax);
            end
                
            
            %center will be at min of TI
            [xpks,xlocs,w,proms]=findpeaks(-TIx);
            %figure;findpeaks(-TIx,xshifts,'Annotate','extents')
            %title('x center')
            [~,nxc]=max(proms);
            bxc=xshifts(xlocs(nxc)); %xshift in pixels
            
            [ypks,ylocs,w,proms]=findpeaks(-TIy);
            %figure;findpeaks(-TIy,yshifts,'Annotate','extents')
            %title('y center')
            [~,nyc]=max(proms);
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
            Atb=[1,0,0;0,1,0;b(1),-b(2),1];
            AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
            
            Arx=[-1,0,0;0,1,0;0,0,1];
            Ary=[1,0,0;0,-1,0;0,0,1];
            
            Nthmax=361;
            thetavals=linspace(0,360,Nthmax)*pi/180;
            
            
            TIrotr=zeros(Nthmax,1);
            TIrefr=zeros(Nthmax,1);
            
            for llth=1:Nthmax
            
                theta=thetavals(llth);
                Arot=[cos(theta), -sin(theta), 0; sin(theta), cos(theta), 0; 0, 0, 1];
                ArotI=[cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
            
                
                %rotate by theta
                tform=AtbI*Arot*Atb;
                TIrotr(llth)=transinfo(imdata,tform,Dthresh,Pmax);
                
                %reflect about axis at theta
                tform=AtbI*ArotI*Ary*Arot*Atb; %bring axis to x-axis, reflect y->-y, rotate back
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
            Atb=[1,0,0;0,1,0;b(1),-b(2),1];
            AtbI=[1,0,0;0,1,0;-b(1),b(2),1];
            
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
                tform=AtbI*Arot*Atb;
                TIrotb(llth)=transinfo(imdata,tform,Dthresh,Pmax);
                
                %reflect about axis at theta
                tform=AtbI*ArotI*Ary*Arot*Atb; %bring axis to x-axis, reflect y->-y, rotate back
                TIrefb(llth)=transinfo(imdata,tform,Dthresh,Pmax);
                
               
                
                 
                
             end
                
            % figure;
            % subplot(1,2,1),polarplot(thetavals,TIrotb'),title('rotation');
            % subplot(1,2,2),polarplot(thetavals,TIrefb'),title('reflection');
            
            xClim=pxc+xshifts([1 end]);
            yClim=pyc-yshifts([end 1]);
            yClimR=yshifts([1 end]);
            %ax=subplot(1,3,1);
            
            
            h=figure('Visible', 'off'); %creating figure without displaying it
            subplot(3,2,1), imshow(imdata0);
            %subplot(3,2,2), imagesc(xcoord,ycoord,imdata), colormap(gray), axis(equal), hold on, plot(rxc,ryc,'ro'), plot(bxc,byc,'c*'), hold off;
            ax=subplot(3,2,2); hold on, image(imdata), colormap(gray); axis image, axis off; 
             set(ax,'YDir','reverse','DataAspectRatio',[1 1 1]);
             xlim(xClim); ylim(yClim);
            plot((pxc+rxc),(pyc-ryc),'ro');
            plot((pxc+bxc),(pyc-byc),'c*');
            hold off
            subplot(3,2,3),plot(thetavals*180/pi,TIrotr'),title('rotation - ctr rot (red)');
            subplot(3,2,4),plot(thetavals*180/pi,TIrefr'),title('reflection - ctr rot (red)');
            subplot(3,2,5),plot(thetavals*180/pi,TIrotb'),title('rotation - ctr ref (cyan)');
            subplot(3,2,6),plot(thetavals*180/pi,TIrefb'),title('reflection - ctr ref (cyan)');
            %print('FlowerSymmetries.png','-dpng')

            %%Writing to folder

            dataFolder='data';
            plotFolder='plots';

            if ~exist(dataFolder,'dir')
                mkdir(dataFolder)
            end

            if ~exist(plotFolder,'dir')
                mkdir(plotFolder)
            end
            
            %pwd for "path to working directory"
            topPath=fullfile(pwd,topFolder);
            
            %Gets path removing topFolder/ (e.g. angiosperms/)
            relativeFolder=erase(file.folder, [topPath,filesep]);

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