source("jngcaFunctions.R")
### Outcome

load("EstimatedComponents_example.Rda")

# Create all the plots for joint components
out_true1 = plotNetwork(Sytrue[,1], title='Truth',qmin=0.005, qmax=0.995,path = 'community_affiliation_mmpplus_example.csv') 

out_true2 = plotNetwork(Sytrue[,2], title='Truth',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus_example.csv') 

out_joint1 = plotNetwork(SyjointICA[,1], title='Joint ICA',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus_example.csv') 

out_joint2 = plotNetwork(SyjointICA[,2], title='Joint ICA',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus_example.csv') 

out_rhoLarge1 = plotNetwork(Sy_rhoLarge[,1], title='large~rho',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus_example.csv') 

out_rhoLarge2 = plotNetwork(Sy_rhoLarge[,2], title='large~rho',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus_example.csv')

out_mcca1 = plotNetwork(SymCCA[,1], title='mCCA+jICA',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus_example.csv') 

out_mcca2 = plotNetwork(SymCCA[,2], title='mCCA+jICA',qmin=0.005, qmax=0.995, path = 'community_affiliation_mmpplus_example.csv') 
#### there are errors in this step
##### Stop at here

par(mfrow = c(2,4))
lgrid=50
image(vec2net(Sytrue[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n") #Truth
image(vec2net(Sytrue[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n") #Truth
image(vec2net(Sy_rhoLarge[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n") # large rho
image(vec2net(Sy_rhoLarge[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n") # large rho
image(vec2net(SyjointICA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n") # Joint ICA
image(vec2net(SyjointICA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n") #Joint ICA
image(vec2net(SymCCA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n") #mCCA+jICA
image(vec2net(SymCCA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n") #mCCA+jICA



# Compare loadings from the components of Y shared structure
############################################################################
mmp_modules = read.csv('community_affiliation_mmpplus.csv')
mmp_order = order(mmp_modules$Community_Vector)
Community = factor(mmp_modules$Community_Label)[mmp_order]

# 1st component loadings
dataLoad <- data.frame(x = rep(c(1:379), 4), loadingsum2 = c(out_true1$loadingsummary[mmp_order], out_rhoLarge1$loadingsummary[mmp_order], out_joint1$loadingsummary[mmp_order], out_mcca1$loadingsummary[mmp_order]), Community = rep(Community, 4), Method = rep(c("Truth", "SING", "Joint~ICA", "mCCA+jICA"), each = 379))

dataLoad$Method <- relevel(as.factor(dataLoad$Method), "Truth")

# 2nd component loadings
dataLoad2 <- data.frame(x = rep(c(1:379), 4), loadingsum2 = c(out_true2$loadingsummary[mmp_order], out_rhoLarge2$loadingsummary[mmp_order], out_joint2$loadingsummary[mmp_order], out_mcca2$loadingsummary[mmp_order]), Community = rep(Community, 4), Method = rep(c("Truth", "SING", "Joint~ICA",  "mCCA+jICA"), each = 379))

dataLoad2$Method <- relevel(as.factor(dataLoad2$Method), "Truth")

dataLoad$component <- "Component~1"
dataLoad2$component <- "Component~2"

dataAll <- rbind(dataLoad, dataLoad2)

dataAll$Method= factor(dataAll$Method, levels = c("Truth", "SING", "Joint~ICA",  "mCCA+jICA"))

# Relevle

pdf(file = "LargeScaleBothComponentsY_mCCA.pdf", width = 14, height = 4)
p = ggplot(dataAll, aes(x = x, y = loadingsum2, col = Community))+geom_point(size = 4) + facet_grid(component~Method, labeller = label_parsed)+xlab('MMP Index')+ylab('L1 Norm of the Rows') + theme(text = element_text(size=20))
print(p)
dev.off()


# Compare visual networks from the components of Y shared structure
############################################################################
qmin=0.005
qmax=0.995
labels = c('VI','SM','DS','VS','DM','CE','SC')
coords = c(0,70.5,124.5,148.5,197.5,293.5,360.5)
zmin = -6
zmax = 6

meltsub = create.graph.long(out_true1$netmat,mmp_order)
meltsub$Method = "Truth"
meltsub$component = "Component~1"
meltsubT = create.graph.long(out_true2$netmat,mmp_order)
meltsubT$Method = "Truth"
meltsubT$component = "Component~2"
meltsub = rbind(meltsub, meltsubT)


meltsub2 = create.graph.long(out_joint1$netmat, mmp_order)
meltsub2$Method = "Joint~ICA"
meltsub2$component = "Component~1"
meltsub2T = create.graph.long(out_joint2$netmat,mmp_order)
meltsub2T$Method = "Joint~ICA"
meltsub2T$component = "Component~2"
meltsub2 = rbind(meltsub2, meltsub2T)


meltsub3 = create.graph.long(out_rhoLarge1$netmat, mmp_order)
meltsub3$Method = "Large~rho"
meltsub3$component = "Component~1"
meltsub3T = create.graph.long(out_rhoLarge2$netmat,mmp_order)
meltsub3T$Method = "Large~rho"
meltsub3T$component = "Component~2"
meltsub3 = rbind(meltsub3, meltsub3T)

meltsub4 = create.graph.long(out_mcca1$netmat, mmp_order)
meltsub4$Method = "mCCA+jICA"
meltsub4$component = "Component~1"
meltsub4T = create.graph.long(out_mcca2$netmat,mmp_order)
meltsub4T$Method = "mCCA+jICA"
meltsub4T$component = "Component~2"
meltsub4 = rbind(meltsub4, meltsub4T)


meltsubAll = rbind(meltsub, meltsub2, meltsub4, meltsub3)
meltsubAll$Method <- as.factor(meltsubAll$Method)
meltsubAll$Method <- relevel(meltsubAll$Method, "Truth")

meltsubAll$Method= factor(meltsubAll$Method, levels = c("Truth", "Joint~ICA",  "mCCA+jICA", "Large~rho"))

# Truncate small values for ease of signal visualization
meltsubAll$value[abs(meltsubAll$value) < 1] = 0

g2 = ggplot(meltsubAll, aes(X1, X2,fill=value)) + geom_tile()+ scale_fill_gradient2(low = "blue",  high = "red", limits=c(zmin,zmax), oob=squish)+labs(x = "Node 1", y = "Node 2") + coord_cartesian(clip='off',xlim=c(-0,390)) + facet_grid(component~Method, labeller = label_parsed)

for (i in 1:7) {
  if (i!=3) {
    g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+10), xmax = (coords[i]+10), ymin = -7, ymax = -7)
  } else{
    g2 = g2+geom_hline(yintercept = coords[i],linetype="dotted",size=0.5)+geom_vline(xintercept = coords[i],linetype="dotted",size=0.5)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),ymin = (coords[i]+10), ymax = (coords[i]+10), xmin = 385, xmax = 385)+annotation_custom(grob = textGrob(label = labels[i], hjust = 0, gp = gpar(cex = 1)),xmin = (coords[i]+1), xmax = (coords[i]+1), ymin = -7, ymax = -7)
  }
}

g2  = g2 + theme(text = element_text(size=24))

g2  = g2 + theme(text = element_text(size=30))
print(g2)

g2  = g2 + theme(text = element_text(size=40))
print(g2)

# pdf(file = "LargeScaleComponentsYnetwork.pdf", width = 14, height = 6)
# print(g2)
# dev.off()

png(filename = "LargeScaleComponentsYnetwork_mCCA.png",
    width = 2500, height = 2000)
print(g2)
dev.off()


# Create cifti files:
#############################################################################################
source('makecifti.R')

# match joint ICA:
SxjointICA_matched = matchICA(S = SxjointICA,template = Sxtrue)
cor(SxjointICA_matched,Sxtrue)
# joint components do not contain the joint signal
cor(Sx_rhoSmall,Sxtrue)

apply(Sx_rhoSmall,2,function(x) mean(x^3))
apply(SxjointICA_matched,2,function(x) mean(x^3))

# If necessary run the following:
#Sxtrue = signchange(Sxtrue)
#Sx_rhoLarge = signchange(Sx_rhoLarge)
#SxjointICA = signchange(SxjointICA)
#apply(t(sIxs),2,function(x) mean(x^3))
#SIxs = signchange(t(sIxs))
# This function has dependencies that will need to be modified. It requires matlab and wb_command:
makecifti(Sxtrue,"Sxtrue.dtseries.nii")
makecifti(Sx_rhoSmall,"Sx_rhoSmall.dtseries.nii")
makecifti(SxjointICA_matched,"SxjointICA.dtseries.nii")

# Alternative showing how to specify pathnames on other systems:
#makecifti(Sxtrue,"Sxtrue.dtseries.nii", wbloc='/home/benjamin/Applications2/workbench/bin_linux64/wb_command', matlabloc='/usr/local/bin/matlab')

