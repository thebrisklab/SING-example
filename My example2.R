source("jngcaFunctions.R")
### Outcome

load("EstimatedComponents_example.Rda")

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






par(mfrow = c(2,4))
image(vec2net(Sytrue[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sytrue[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sy_rhoSmall[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Large") # large rho
image(vec2net(Sy_rhoSmall[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Large") # large rho
image(vec2net(SyjointICA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA") # Joint ICA
image(vec2net(SyjointICA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA") #Joint ICA
image(vec2net(SymCCA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA") #mCCA+jICA
image(vec2net(SymCCA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA") #mCCA+jICA








