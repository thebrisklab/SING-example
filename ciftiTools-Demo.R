### Load `ciftiTools`
# install.packages("ciftiTools")
library(ciftiTools)

### Load the Connectome Workbench.
ciftiTools.setOption("wb_path", "C:/Software/workbench")


### Make example matrix data.
nVertLeft <- 5762
nVertRight <- 5762
nTime <- 400

dataLeft <- matrix(rnorm(nVertLeft*nTime), ncol=nTime)
dataRight <- matrix(rnorm(nVertRight*nTime), ncol=nTime)

### Construct a `"xifti"` object.
xii <- as.xifti(cortexL=dataLeft, cortexR=dataRight)

### Optional: plot the first column.
# plot(xii, title="First column")

### Write out the `"xifti"` object.
write_cifti(xii, "my_cifti.dtseries.nii")
view_xifti_surface(xii)




dataLeft <- Sxtrue[1:29696,1:2]
dataRight <- Sxtrue[29697:59412,1:2]
xii <- as.xifti(cortexL=dataLeft, cortexR=dataRight,s)
a=as.matrix.xifti(dataLeft)


matrix(rep(0,728*2),728,2)

