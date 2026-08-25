%Scratch work to show that adding noise doesn't change TI curve much

clear

%Adding path for auxiliary functions
parent=fileparts(fileparts(pwd));
addpath(fullfile(parent,'AuxFunctions'));

%Read and store image
tmp=imread("figWhiteBack.jpg");
imdata0=tmp;

imdata=(255-im2gray(tmp)); %Convert to grayscale, take negative

%Add random Gaussian noise
imdata=imnoise(imdata,'gaussian',0,0.01); %Converts to doubles on [0,1] plus noise
%imdata=imnoise(imdata,'salt & pepper',0.05);
%imdata=im2uint8(imdata)+1;
%imdata=imgaussfilt(im2uint8(imdata)+1,1);


[M,N]=size(imdata);

Pmax=256;
Dthresh=15;

% center of image 
pyc=ceil(M/2); 
pxc=ceil(N/2); 

%index vectors
%(1,1) is at the top left of the image 
% coordinates to search through rectangle in the middle of the image
xcoord=(pxc-floor(M/8)):(pxc+ceil(M/8)); %positive right, same as indices
ycoord=(pyc-floor(N/8)):(pyc+ceil(N/8)); %positive up, backwards from indices


%figure;imshow(imdata)

%%%%%%%%%%%%%
%Find Center
%%%%%%%%%%%%%

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
            TIcenter(llx,kky,mth)=transinfo(imdata(1:npskip:end,1:npskip:end),Atform,Dthresh,Pmax);
        end

    end

end

%figure;imagesc(xshifts,yshifts,TIcenter(:,:,1)')
%figure;imagesc(xshifts,yshifts,(max(TIcenter,[],3)-min(TIcenter,[],3))')
%figure;imagesc(xshifts,yshifts,sum(abs(diff(TIcenter,[],3)).^2,3)')

%arr=sum(abs(diff(TIcenter,[],3)).^2,3)';

%For each center point in the grid, find the minimum
for rarr=1:NYmax
    for carr=1:NXmax
        arr(rarr,carr)=min(TIcenter(rarr,carr,:));
    end
end

%arr=TIcenter;
minimum=min(min(arr)); %Find the minimum of the mins
[xidmax,yidmax]=find(arr==minimum,1); %index of center
%xidmax=mod(xidmaxtmp,NXmax);
rxc=xshifts(xidmax); %center value in relation to midpoint
ryc=yshifts(yidmax);
%rxc=xshifts(11);
%ryc=yshifts(13);

%b=[bxc,byc];
b=[rxc,ryc];
%b=[xFc,yFc]; %geometric centre of flower



%% Compute TI curves about center

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


figure;
plot(thetavals,TIrotr,'LineWidth',3.0);
xlim([0,2*pi]);

figure;
plot(thetavals,TIrefr,'Linewidth',3.0);
xlim([0,2*pi]);

beep