# TransInfo
Matlab script to automatically compute the symmetry of image data using *transformation information (TI)*. The symmetries are then compared using either principle component analysis (PCA) or a computationally tractable version of Gromov-Wasserstein distance, known as GW-tau.

## User instructions:
1. Place your images in a folder with some label. We have labeled our image file "angiosperms"; however, any name such as "Images" or "ExperimentalData" is okay. Just make sure to avoid the name of an existing file (i.e., "data", "plots", "pcaData", "pcaPlots", and "all_pcaPlots"). 

There may be subfolders as well to classify data that we wish to compare. In our example, we wanted to compare different classes of flowers within the phylum angiosperms, so we placed each class in a different folder. If you want to compare experimental data, you may have subfolders with images from different experimental conditions.

2. Download the code, and place the image folder created in step 1 into the same directory/folder.

3. Run the file with title 'user_input.m'. It will prompt you to input some information such as the name of your image file and...

Relevant data will be output in the following folders.
- all_pcaPlots: This folder contains a plot of the first two principle components of the TI for all of your images. The data will be color coded based on the file that you used for classification.

- pcaPlots: This folder contains plots of the first two principle components of the tfor each subfolder (PCA computed separately for each class). The folder scheme within the pcaPlots folder will match your original folder scheme for the image data.

- pcaData: This folder contains the pcaData for all of your image classes. The PCA is computed separately for each class.

- plots: This folder contains images of the TI curves for each image. The folder organization follows the scheme set in your image data folder.

- data: Contains the raw TI data used for "plots."


## Files:

- angiosperms: Folder of images to analyze

- transifo.m: Function for computing transformation information associated to a given transformation

- TransInfo_images.m: File to compute the center and TI curve of each image

- class_pca.m: Runs the PCA within each class given by the subfolders of the original image folder, "angiosperms."

- ti_pca_all.m: Computes the PCA of all image data, with the scatter plot classifying the images by color.



## References
Thank you to Dr. Punit Gandhi at Virginia Commonwealth University for providing the original code, which was forked. Code and data added to the forked repository was created by Dr. Adriana Dawes at The Ohio State University, with minor modifications done by me.
