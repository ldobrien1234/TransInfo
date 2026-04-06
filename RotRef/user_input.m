%Gets user inputs and completes data analysis

clear

%Selecting data folder from parent directory
topFolder=uigetdir(fileparts(pwd), 'Please select your data folder.');
topFolder=erase(topFolder, [fileparts(pwd),filesep]);


%Computes TI curves for the input folder and classifies the data
%appropriately

TransInfo_images(topFolder);

%Computes the PCA by class, given by the folders within the topFolder
%Plots the TI curves along the first two principle components
class_pca();

%Computes the PCA for all TI curves (all classes) together
%Plots the TI curves along the first two principle components
ti_pca_all();

%Uses GWtau and multi-dimensional scaling to plot TI curves in scatter plot
gwtau_embed_all();

beep %makes sound when finished