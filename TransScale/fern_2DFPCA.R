# FPCA and embedding on 2D fern TI curves

# Installing the required packages
required_packages <- c("dplyr", "fda", "lsei", "raster")
(missed_packages <- setdiff(required_packages, rownames(installed.packages())))
if (length(missed_packages)) {
  sapply(missed_packages, install.packages)
}
library(R.matlab) #to load data from MATLAB file
library(dplyr)
library(fda)
library(lsei)

# Set the working directory to the folder with the files "Mnist.RData" and "funs_2DFPCA.R"
setwd('/Users/liamobrien/Documents/Dawes Lab/TI Project/Flowers/TransInfo2/TransScale')

# Load the data and source code
source("funs_2DFPCA.R") #Runs R code from specified file
mat_data=readMat("fPCA_vars.mat")


observed <- mat_data$TI.frames.pca
timepoints1 <- mat_data$pca.xind
timepoints2 <- mat_data$pca.yind

#Change data types from matrix/array of lists to a list of lists
observed <- lapply(observed, function(cell_item) {
  # as.numeric removes the MATLAB 'matrix' attributes from the vector
  return(as.numeric(unlist(cell_item)))
})
timepoints1<- lapply(timepoints1, function(cell_item) {
  # as.numeric removes the MATLAB 'matrix' attributes from the vector
  return(as.numeric(unlist(cell_item)))
})

timepoints2 <- lapply(timepoints2, function(cell_item) {
  # as.numeric removes the MATLAB 'matrix' attributes from the vector
  return(as.numeric(unlist(cell_item)))
})

# #Points in the mesh of the TI curve
# timepoints1<-rep(list(1:30),length(observed))
# timepoints2<-rep(list(1:30),length(observed))


# Center the data 

#Observed gives images as lists of pixel intensities
#Timepoints 1 and 2 give first and second index for each pixel (hence the 2D structure)
meanobserve <- Reduce("+", observed)/length(observed[[1]])
observed <- lapply(observed, function(x) {
  x-meanobserve
})


# Initialize the spline parameters 
maxpixelsize <- 30
nbasis <- 12
library(fda)
beta1 <- matrix(1, nrow = nbasis, ncol = nbasis)
for (i in 1:nbasis) {
  beta1[, i] <- seq(0.01, 0.12, by = 0.01)
}

basis1 <- create.bspline.basis(rangeval = c(1, maxpixelsize), nbasis = nbasis, norder = 4)
basis2 <- create.bspline.basis(rangeval = c(1, maxpixelsize), nbasis = nbasis, norder = 4)

initializeGlobalXmat(timepoints1, timepoints2, basis1, basis2)
	
previous_beta <- list()
pc_list <- list()
result_list <- list()




# Fit the first 3 FPCs
numFPC <- 3
for (i in 1:numFPC) {
  if (i == 1) {
	
	 # First FPC
    res_first <- first_FPC_2d_image(beta1, observed, timepoints1, timepoints2, basis1, basis2, threshold = 1e-4)
    result_list[[1]] <- res_first
    previous_beta[[1]] <- res_first$beta #Spline coefficient
    pc_list[[1]] <- res_first$pc_fit
  } else {
    
	 # Higher order FPCs
    res_higherorder <- second_FPC_conditional_2d_image(beta1, pc_index = i, observed, timepoints1, timepoints2, basis1, basis2, betalist = previous_beta, threshold = 1e-4)
    result_list[[i]] <- res_higherorder
    previous_beta[[i]] <- res_higherorder$beta
    pc_list[[i]] <- res_higherorder$pc_fit
    scores <- res_higherorder$sfit
  }
}

library(raster)

# # Plot the first 3 FPCs	for images of digits "0" and "1" in the MNIST dataset
# quartz() # open a new plot window (for Mac OS only)
# par(mfrow = c(2,2))
# for (i in 1:numFPC){
# 
#   plotindex <- i
#   res <- result_list[[plotindex]]
#   sfit <- result_list[[length(result_list)]]$sfit
# 
#   plotmat = matrix(0, ncol = 28, nrow = 28)
#   for (x in 1:28){
#     for (y in 1:28){
#       #fd2d gives tensor spline basis function from coefficients beta
#       plotmat[x, y] <- (eval.fd2d(x, y, res$pc_fit) * 1000)
#     }
#   }
# 
#   library(raster)
#   plot(raster(plotmat), col = grey.colors(10, start = 0, end = 1), main = paste0("FPC ", plotindex))
# }

  
  
  
  
# library(ggplot2)
#
# df <- data.frame(
#   fpc1 = scores[,1],
#   fpc2 = scores[,2]
# )
#
#   ggplot(df, aes(x = fpc1, y = fpc2)) +
#     geom_point(alpha = 0.6, size = 2) +
#     theme_minimal() +
#     labs(title = "Functional PCA Scores",
#          x = "First Principal Component (Score)",
#          y = "Second Principal Component (Score)")
  



#Same scatter plot as above formatted to look like MATLAB
library(ggplot2)

matlab_blue <- "#0072BD"

df <- data.frame(
  fpc1 = scores[,1],
  fpc2 = scores[,2]
)

ggplot(df, aes(x = fpc1, y = fpc2)) +
  geom_point(color = matlab_blue, alpha = 1, size = 2.3) +
  labs(title = "Functional PCA Scores",
       x = "First Principal Component (Score)",
       y = "Second Principal Component (Score)") +
  theme_minimal() +
  theme(
    # 1. Remove gridlines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),

    # 2. L-shaped border
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),

    # 3. MOVE TICKS INSIDE (The negative value trick)
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(-0.15, "cm"),

    # 4. PUSH LABELS AWAY (Prevents them from touching the inward ticks)
    axis.text.x = element_text(margin = margin(t = 10)), # t = top margin
    axis.text.y = element_text(margin = margin(r = 10)), # r = right margin

    # 5. Typography
    text = element_text(family = "sans", size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )