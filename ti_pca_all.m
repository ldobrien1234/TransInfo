close all
clear

%Load data file
load('ti_early_angiosperms.mat')
%Rows are TI data for a given flower, columns are theta mesh points
[rws,cols]=size(tirotr_all); %Size of TI data function
class_count(1)=rws; %Number of flowers in 1st class (early angiosperms)
tirotr_angiosperms(:,:)=tirotr_all;
tirefr_angiosperms(:,:)=tirefr_all;
tirotb_angiosperms(:,:)=tirotb_all;
tirefb_angiosperms(:,:)=tirefb_all;


%Append flowers from second data class to bottom of matrices constructed
%above; again, each row is a flower's TI data
%When loading the new file, variables tirotr_all, tirefr_all, tirotb_all,
%tirefb_all are overwritten
load('ti_eudicots_asterids.mat')
[rws,cols]=size(tirotr_all);
class_count(2)=rws;
tirotr_angiosperms(class_count(1)+1:sum(class_count(1:2)),:)=tirotr_all;
tirefr_angiosperms(class_count(1)+1:sum(class_count(1:2)),:)=tirefr_all;
tirotb_angiosperms(class_count(1)+1:sum(class_count(1:2)),:)=tirotb_all;
tirefb_angiosperms(class_count(1)+1:sum(class_count(1:2)),:)=tirefb_all;

load('ti_eudicots_early_diverging_eudicots.mat')
[rws,cols]=size(tirotr_all);
class_count(3)=rws;
tirotr_angiosperms(sum(class_count(1:2))+1:sum(class_count(1:3)),:)=tirotr_all;
tirefr_angiosperms(sum(class_count(1:2))+1:sum(class_count(1:3)),:)=tirefr_all;
tirotb_angiosperms(sum(class_count(1:2))+1:sum(class_count(1:3)),:)=tirotb_all;
tirefb_angiosperms(sum(class_count(1:2))+1:sum(class_count(1:3)),:)=tirefb_all;

load('ti_eudicots_rosids.mat')
[rws,cols]=size(tirotr_all);
class_count(4)=rws;
tirotr_angiosperms(sum(class_count(1:3))+1:sum(class_count(1:4)),:)=tirotr_all;
tirefr_angiosperms(sum(class_count(1:3))+1:sum(class_count(1:4)),:)=tirefr_all;
tirotb_angiosperms(sum(class_count(1:3))+1:sum(class_count(1:4)),:)=tirotb_all;
tirefb_angiosperms(sum(class_count(1:3))+1:sum(class_count(1:4)),:)=tirefb_all;

load('ti_eudicots_superasterids.mat')
[rws,cols]=size(tirotr_all);
class_count(5)=rws;
tirotr_angiosperms(sum(class_count(1:4))+1:sum(class_count(1:5)),:)=tirotr_all;
tirefr_angiosperms(sum(class_count(1:4))+1:sum(class_count(1:5)),:)=tirefr_all;
tirotb_angiosperms(sum(class_count(1:4))+1:sum(class_count(1:5)),:)=tirotb_all;
tirefb_angiosperms(sum(class_count(1:4))+1:sum(class_count(1:5)),:)=tirefb_all;

load('ti_eudicots_superrosids.mat')
[rws,cols]=size(tirotr_all);
class_count(6)=rws;
tirotr_angiosperms(sum(class_count(1:5))+1:sum(class_count(1:6)),:)=tirotr_all;
tirefr_angiosperms(sum(class_count(1:5))+1:sum(class_count(1:6)),:)=tirefr_all;
tirotb_angiosperms(sum(class_count(1:5))+1:sum(class_count(1:6)),:)=tirotb_all;
tirefb_angiosperms(sum(class_count(1:5))+1:sum(class_count(1:6)),:)=tirefb_all;

load('ti_monocots.mat')
[rws,cols]=size(tirotr_all);
class_count(7)=rws;
tirotr_angiosperms(sum(class_count(1:6))+1:sum(class_count(1:7)),:)=tirotr_all;
tirefr_angiosperms(sum(class_count(1:6))+1:sum(class_count(1:7)),:)=tirefr_all;
tirotb_angiosperms(sum(class_count(1:6))+1:sum(class_count(1:7)),:)=tirotb_all;
tirefb_angiosperms(sum(class_count(1:6))+1:sum(class_count(1:7)),:)=tirefb_all;

sz=50; 

 [coeff,score,latent,tsquared,explained] = pca(tirotr_angiosperms);
h=figure; 
scatter(score(1:class_count(1),1),score(1:class_count(1),2),sz,'red','filled')
hold on
%for i=1:file_count
%   scatter(score(sum(class_count(1:i-1))+1:class_count(i),1),score(sum(class_count(1:i-1))+1:class_count(i),2),sz,'red','filled')
%end
scatter(score(class_count(1)+1:sum(class_count(1:2)),1),score(class_count(1)+1:sum(class_count(1:2)),2),sz,'green','filled')
scatter(score(sum(class_count(1:2))+1:sum(class_count(1:3)),1),score(sum(class_count(1:2))+1:sum(class_count(1:3)),2),sz,'blue','filled')
scatter(score(sum(class_count(1:3))+1:sum(class_count(1:4)),1),score(sum(class_count(1:3))+1:sum(class_count(1:4)),2),sz,'cyan','filled')
scatter(score(sum(class_count(1:4))+1:sum(class_count(1:5)),1),score(sum(class_count(1:4))+1:sum(class_count(1:5)),2),sz,'magenta','filled')
scatter(score(sum(class_count(1:5))+1:sum(class_count(1:6)),1),score(sum(class_count(1:5))+1:sum(class_count(1:6)),2),sz,'yellow','filled')
scatter(score(sum(class_count(1:6))+1:sum(class_count(1:7)),1),score(sum(class_count(1:6))+1:sum(class_count(1:7)),2),sz,'black','filled')

title({
    ['TI by rotation, Rot ctr' ]
    [num2str(explained(1)+explained(2)) ' variance explained' ]
    });
print(h,'-djpeg','tirotr_all_angiosperms_pca.jpg')

 [coeff,score,latent,tsquared,explained] = pca(tirefr_angiosperms);
h=figure;
scatter(score(1:class_count(1),1),score(1:class_count(1),2),sz,'red','filled')
hold on
%for i=1:file_count
%   scatter(score(sum(class_count(1:i-1))+1:class_count(i),1),score(sum(class_count(1:i-1))+1:class_count(i),2),sz,'red','filled')
%end
scatter(score(class_count(1)+1:sum(class_count(1:2)),1),score(class_count(1)+1:sum(class_count(1:2)),2),sz,'green','filled')
scatter(score(sum(class_count(1:2))+1:sum(class_count(1:3)),1),score(sum(class_count(1:2))+1:sum(class_count(1:3)),2),sz,'blue','filled')
scatter(score(sum(class_count(1:3))+1:sum(class_count(1:4)),1),score(sum(class_count(1:3))+1:sum(class_count(1:4)),2),sz,'cyan','filled')
scatter(score(sum(class_count(1:4))+1:sum(class_count(1:5)),1),score(sum(class_count(1:4))+1:sum(class_count(1:5)),2),sz,'magenta','filled')
scatter(score(sum(class_count(1:5))+1:sum(class_count(1:6)),1),score(sum(class_count(1:5))+1:sum(class_count(1:6)),2),sz,'yellow','filled')
scatter(score(sum(class_count(1:6))+1:sum(class_count(1:7)),1),score(sum(class_count(1:6))+1:sum(class_count(1:7)),2),sz,'black','filled')

title({
    ['TI by reflection, Rot ctr' ]
    [num2str(explained(1)+explained(2)) ' variance explained' ]
    });
print(h,'-djpeg','tirefr_all_angiosperms_pca.jpg')

 [coeff,score,latent,tsquared,explained] = pca(tirotb_angiosperms);
h=figure;
scatter(score(1:class_count(1),1),score(1:class_count(1),2),sz,'red','filled')
hold on
%for i=1:file_count
%   scatter(score(sum(class_count(1:i-1))+1:class_count(i),1),score(sum(class_count(1:i-1))+1:class_count(i),2),sz,'red','filled')
%end
scatter(score(class_count(1)+1:sum(class_count(1:2)),1),score(class_count(1)+1:sum(class_count(1:2)),2),sz,'green','filled')
scatter(score(sum(class_count(1:2))+1:sum(class_count(1:3)),1),score(sum(class_count(1:2))+1:sum(class_count(1:3)),2),sz,'blue','filled')
scatter(score(sum(class_count(1:3))+1:sum(class_count(1:4)),1),score(sum(class_count(1:3))+1:sum(class_count(1:4)),2),sz,'cyan','filled')
scatter(score(sum(class_count(1:4))+1:sum(class_count(1:5)),1),score(sum(class_count(1:4))+1:sum(class_count(1:5)),2),sz,'magenta','filled')
scatter(score(sum(class_count(1:5))+1:sum(class_count(1:6)),1),score(sum(class_count(1:5))+1:sum(class_count(1:6)),2),sz,'yellow','filled')
scatter(score(sum(class_count(1:6))+1:sum(class_count(1:7)),1),score(sum(class_count(1:6))+1:sum(class_count(1:7)),2),sz,'black','filled')

title({
    ['TI by rotation, Ref ctr' ]
    [num2str(explained(1)+explained(2)) ' variance explained' ]
    });
print(h,'-djpeg','tirotb_all_angiosperms_pca.jpg')

 [coeff,score,latent,tsquared,explained] = pca(tirefb_angiosperms);
h=figure;
scatter(score(1:class_count(1),1),score(1:class_count(1),2),sz,'red','filled')
hold on
%for i=1:file_count
%   scatter(score(sum(class_count(1:i-1))+1:class_count(i),1),score(sum(class_count(1:i-1))+1:class_count(i),2),sz,'red','filled')
%end
scatter(score(class_count(1)+1:sum(class_count(1:2)),1),score(class_count(1)+1:sum(class_count(1:2)),2),sz,'green','filled')
scatter(score(sum(class_count(1:2))+1:sum(class_count(1:3)),1),score(sum(class_count(1:2))+1:sum(class_count(1:3)),2),sz,'blue','filled')
scatter(score(sum(class_count(1:3))+1:sum(class_count(1:4)),1),score(sum(class_count(1:3))+1:sum(class_count(1:4)),2),sz,'cyan','filled')
scatter(score(sum(class_count(1:4))+1:sum(class_count(1:5)),1),score(sum(class_count(1:4))+1:sum(class_count(1:5)),2),sz,'magenta','filled')
scatter(score(sum(class_count(1:5))+1:sum(class_count(1:6)),1),score(sum(class_count(1:5))+1:sum(class_count(1:6)),2),sz,'yellow','filled')
scatter(score(sum(class_count(1:6))+1:sum(class_count(1:7)),1),score(sum(class_count(1:6))+1:sum(class_count(1:7)),2),sz,'black','filled')

% title({
%     ['TI by reflection, Ref ctr' ]
%     [num2str(explained(1)+explained(2)) ' variance explained' ]
%     });
print(h,'-djpeg','tirefb_all_angiosperms_pca.jpg')













