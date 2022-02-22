### Quick start

# Load the package and point to the Connectome Workbench --------
library(ciftiTools)
ciftiTools.setOption("wb_path", "C:/Software/workbench")

# Read and visualize a CIFTI file -------------------------------
cifti_fname <- ciftiTools::ciftiTools.files()$cifti["dtseries"]
surfL_fname <- ciftiTools.files()$surf["left"]
surfR_fname <- ciftiTools.files()$surf["right"]


xii <- read_cifti(
  cifti_fname, brainstructures="all", 
  surfL_fname=surfL_fname, surfR_fname=surfR_fname,
  resamp_res=4000
)

view_xifti_surface(xii) # or plot(xii)
#view_xifti_volume(xii) if subcortex is present

# Access CIFTI data ---------------------------------------------
cortexL <- xii$data$cortex_left
cortexL_mwall <- xii$meta$medial_wall_mask$left
cortexR <- xii$data$cortex_right
cortexR_mwall <- xii$meta$medial_wall_mask$right
# subcortVol <- xii$data$subcort
# subcortLabs <- xii$meta$subcort$labels
# subcortMask <- xii$meta$subcort$mask
surfL <- xii$surf$cortex_left
surfR <- xii$surf$cortex_right

# Create a `"xifti"` from data ----------------------------------
xii2 <- as.xifti(
  cortexL=cortexL, cortexL_mwall=cortexL_mwall,
  cortexR=cortexR, cortexR_mwall=cortexR_mwall,
  #subcortVol=subcortVol, subcortLabs=subcortLabs,
  #subcortMask=subcortMask,
  #surfL=surfL, surfR=surfR
)

# Write a CIFTI file --------------------------------------------
write_cifti(xii2, "my_cifti.dtseries.nii")


### cifti demo
library(ciftiTools)
ciftiTools.setOption("wb_path", "C:/Software/workbench")

cifti_fnames <- ciftiTools.files()$cifti
surfL_fname <- ciftiTools.files()$surf["left"]
surfR_fname <- ciftiTools.files()$surf["right"]

basename(cifti_fnames["dtseries"])

xii <- read_xifti(cifti_fnames["dtseries"])
xii # same as `summary(xii)`

xii2 <- read_xifti(cifti_fnames["dtseries"], surfL_fname=surfL_fname, surfR_fname=surfR_fname)
all.equal(xii, xii2) # same result

xii_info <- ciftiTools::info_cifti(cifti_fnames["dscalar"])
str(xii_info, nchar.max=50) # shows header structure

read_xifti(cifti_fnames["dtseries"], idx=2) # second column only


out_dir <- "output"

write_xifti(
  xii, 
  file.path(out_dir, "my_cifti.dtseries.nii"), 
  file.path(out_dir, "my_L.surf.gii"), 
  file.path(out_dir, "my_R.surf.gii")
)

# Use default names for everything except left cortex
separated_fnames = separate_cifti(
  cifti_fnames["dscalar_ones"], brainstructures="all", 
  cortexL_fname="my_left_cortex.func.gii", write_dir = out_dir
)


# Files written to `out_dir`, or current working dir. if not specified
basename(separated_fnames)


library(rgl)
rgl::setupKnitr()

# Sometimes the first OpenGL window does not render properly.
rgl::rgl.open()
rgl::rgl.close()

# Surface Visualization
view_xifti_surface(xii)


# Normally `cex.title` doesn't need to be set, as it defaults to a good choice.
#   But when knitting static images this way, the default becomes a bit too big
#   based on how knitting works.
view_xifti_surface(xii, idx=1, zlim=c(1,2), title='color_mode = "sequential"', cex.title=1.3)



xii <- read_cifti(cifti_fnames["dscalar"]) # no GIFTI included, so the default inflated surface is used.
view_xifti_surface(
  xii, idx=1:2, zlim=c(0,5), color_mode = "diverging",
  title='color_mode = "diverging"', cex.title=1.3
)



view_xifti_surface(
  read_cifti(cifti_fnames["dlabel"]), 
  # Interactively, a color legend that displays the label names will also be printed.
  legend_ncol=5, 
  title='color_mode = "qualitative"', cex.title=1.3
)

xiiL <- read_xifti(cifti_fnames["dscalar"], brainstructures="left")
dim(as.matrix(xiiL))


gmeans <- apply_xifti(xiiL, 1, mean)
gmeans


cquants <- apply_xifti(xiiL, 2, quantile, c(.1, .5)) # Quantiles of each column
cquants


xiiR <- read_xifti(cifti_fnames["dscalar"], brainstructures="right")
xii <- combine_xifti(xiiL, xiiR)
xii

convert_xifti(xii, "dtseries") # Convert from dscalar to dtseries (don't use it)

xii <- merge_xifti(xii, xii)
xii

xii <- remove_xifti(xii, "cortex_left") # Now only the right cortex data is included.
xii


xii$meta$cifti$names <- paste("Column", seq(4))
xii <- select_xifti(xii, c(4,3,2)) # Reverse column order & drop the first.
xii


max(as.matrix(xii))

xii <- 1 - exp(xii) / (xii * 2 + 3)
max(as.matrix(xii))



# Reading
surf <- read_surf(surfL_fname)
surf


# Writing
write_surf_gifti(surf, file.path(out_dir, "my.L.surf"))

# Visualizing
plot(surf)

# Resample a `"surf"` object
surf <- resample_surf(surf, 2000, "left")
# Resample a GIFTI file
resample_gifti(surfL_fname, file.path(out_dir, "my.L.2k.surf.gii"), "left", resamp_res=2000)

xii <- as.xifti(
  surfL = load_surf("left", "very inflated"),
  surfR = load_surf("right", "midthickness")
)
plot(xii, title = "Left very inflated | Right midthickness")


