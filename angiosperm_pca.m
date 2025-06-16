clear
close all

files_all=dir('*data.mat'); %Gets all files in current folder ending in data.mat
%files=dir('15.png');
file_count_all=1;

for file=files_all'
   load(file.name);
   tirotr_all(file_count_all,:)=TIrotr;
   tirefr_all(file_count_all,:)=TIrefr;
   tirotb_all(file_count_all,:)=TIrotb;
   tirefb_all(file_count_all,:)=TIrefb;
   file_count_all=file_count_all+1;
end

%c = linspace(1,10,file_count-1);
c = jet(file_count_all-1);
sz=50;

 [coeff,score,latent,tsquared,explained] = pca(tirotr_all);
h=figure; scatter(score(:,1),score(:,2),sz,c,'filled')
title({
    ['TI by rotation, Rot ctr' ] 
    [num2str(explained(1)+explained(2)) ' variance explained' ] 
    });
%ftitle=[file.name,'tirotr_all_pca.jpg'];
print(h,'-djpeg','tirotr_all_pca.jpg')

 [coeff,score,latent,tsquared,explained] = pca(tirefr_all);
h=figure; scatter(score(:,1),score(:,2),sz,c,'filled')
title({
    ['TI by reflection, Rot ctr' ] 
    [num2str(explained(1)+explained(2)) ' variance explained' ]
    });
%ftitle=[file.name,'tirotr_all_pca.jpg'];
print(h,'-djpeg','tirefr_all_pca.jpg')


 [coeff,score,latent,tsquared,explained] = pca(tirotb_all);
h=figure; scatter(score(:,1),score(:,2),sz,c,'filled')
title({
    ['TI by rotation, Ref ctr' ] 
    [num2str(explained(1)+explained(2)) ' variance explained' ]
    });
%ftitle=[file.name,'tirotr_all_pca.jpg'];
print(h,'-djpeg','tirotb_all_pca.jpg')


 [coeff,score,latent,tsquared,explained] = pca(tirefb_all);
h=figure; scatter(score(:,1),score(:,2),sz,c,'filled')
title({
    ['TI by reflection, Ref ctr' ] 
    [num2str(explained(1)+explained(2)) ' variance explained' ]
    });
%ftitle=[file.name,'tirotr_all_pca.jpg'];
print(h,'-djpeg','tirefb_all_pca.jpg')
close all

tirotr_early_angiosperms=tirotr_all;
tirefr_early_angiosperms=tirefr_all;
tirotb_early_angiosperms=tirotb_all;
tirefb_early_angiosperms=tirefb_all;

save('ti_early_angiosperms')

