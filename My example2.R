source("jngcaFunctions.R")
### Outcome

load("EstimatedComponents_example.Rda")
load("Estimated joint subject score_example.Rda")

# Create all the plots for joint components
 
#### there are errors in this step
##### Stop at here


lgrid = 33
par(mfrow = c(2,4))

image(matrix(Sxtrue[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth")
image(matrix(Sxtrue[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth")
image(matrix(Sx_rhoSmall[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Small")
image(matrix(Sx_rhoSmall[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Small")
image(matrix(Sx_rhoLarge[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Large")
image(matrix(Sx_rhoLarge[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Large")
image(matrix(Sx_rho0[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=0")
image(matrix(Sx_rho0[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=0")


lgrid = 33
par(mfrow = c(2,4))
image(matrix(Sxtrue[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth")
image(matrix(Sxtrue[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth")
image(matrix(Sx_rhoLarge[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Large")
image(matrix(Sx_rhoLarge[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Large")
image(matrix(SxjointICA[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA")
image(matrix(SxjointICA[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA")
image(matrix(SxmCCA[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA")
image(matrix(SxmCCA[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA")



lgrid = 33
par(mfrow = c(2,4))
image(matrix(Sxtrue[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth")
image(matrix(Sxtrue[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth")
image(matrix(Sx_rhoLarge[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Large")
image(matrix(Sx_rhoLarge[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Large")
image(matrix(SxjointICA[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA")
image(matrix(SxjointICA[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA")
image(matrix(SxmCCA[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA")
image(matrix(SxmCCA[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA")




par(mfrow = c(2,4))
image(vec2net(Sytrue[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sytrue[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sy_rhoSmall[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Large") # large rho
image(vec2net(Sy_rhoSmall[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Large") # large rho
image(vec2net(SyjointICA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA") # Joint ICA
image(vec2net(SyjointICA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA") #Joint ICA
image(vec2net(SymCCA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA") #mCCA+jICA
image(vec2net(SymCCA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA") #mCCA+jICA




par(mfrow = c(2,4))
image(vec2net(Sytrue[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sytrue[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sy_rhoSmall[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Large") # large rho
image(vec2net(Sy_rhoSmall[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Large") # large rho
image(vec2net(SyjointICA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA") # Joint ICA
image(vec2net(SyjointICA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA") #Joint ICA
image(vec2net(SymCCA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA") #mCCA+jICA
image(vec2net(SymCCA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA") #mCCA+jICA


###Figure for the joint subject score

#True mj
ggplot(data = trueMj)+
  geom_point(mapping = aes(y=mj1,x=number))+
  ggtitle("TrueMj,Comp1")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggplot(data = trueMj)+
  geom_point(mapping = aes(y=mj2,x=number))+
  ggtitle("TrueMj,Comp2")+
  theme_bw()+
  theme(panel.grid = element_blank())

#SING mj
ggplot(data = SINGMj)+
  geom_point(mapping = aes(y=mj1,x=number))+
  ggtitle("SINGMj,Comp1")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggplot(data = SINGMj)+
  geom_point(mapping = aes(y=mj2,x=number))+
  ggtitle("SINGmj,Comp2")+
  theme_bw()+
  theme(panel.grid = element_blank())

#ICA mj
ggplot(data = ICAMj)+
  geom_point(mapping = aes(y=mj1,x=number))+
  ggtitle("Joint Score,Comp1")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggplot(data = ICAMj)+
  geom_point(mapping = aes(y=mj2,x=number))+
  ggtitle("Joint Score,Comp2")+
  theme_bw()+
  theme(panel.grid = element_blank())





