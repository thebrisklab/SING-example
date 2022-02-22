### SING example with package
library(SING)

# Generate data, will return both components and mixing matrix
set.seed(0573452)
data <- generateData_v3(nsubject = 48, snr = c(1, 1), vars = c(0.005, 0.005))
save(data,file='data_example.Rda')


lgrid = 33
par(mfrow = c(2,4))
# Plot true joint and individual components for X
image(matrix(data$sjX[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="X joint component1")
image(matrix(data$sjX[2,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="X joint component2")
image(matrix(data$siX[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="X independent component2")
image(matrix(data$siX[2,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="X independent component2")

# Plot true joint and individual components for Y
image(vec2net(data$sjY[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Y joint component1")
image(vec2net(data$sjY[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Y joint component2")
image(vec2net(data$siY[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Y independent component1")
image(vec2net(data$siY[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Y independent component2")

# Data processing steps:
# Center X and Y
n = nrow(data$dX)
pX = ncol(data$dX)
pY = ncol(data$dY)
dXcentered <- data$dX - matrix(rowMeans(data$dX), n, pX, byrow = F)
dYcentered <- data$dY - matrix(rowMeans(data$dY), n, pY, byrow = F)


## Apply separate JB 
#############################################################
# JB on X
estX_JB = lngca(xData = t(data$dX), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
Uxfull <- estX_JB$Ws  # Ax = Ux %*% Lx, where Lx is the whitened matrix from covariance matrix of dX.
Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(data$dX)) # NOTE: for centered X, equivalent to xData %*% sData/(px-1)


# JB on Y
estY_JB = lngca(xData = t(data$dY), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
Uyfull <- estY_JB$Ws
My_JB = est.M.ols(sData = estY_JB$S, xData = t(data$dY))

# Get joint components out
# Greedy Match
matchMxMy = greedymatch(t(Mx_JB), t(My_JB), Ux = t(Uxfull), Uy = t(Uyfull))
permJoint <- permTestJointRank(matchMxMy$Mx,matchMxMy$My) # alpha = 0.01, nperm=1000
pval_joint = permJoint$pvalues
joint_rank = permJoint$rj
joint_rank


# Whiten dX and dY
###################################

# For X
# Scale rowwise
est.sigmaXA = tcrossprod(dXcentered)/(pX-1)  ## dXcentered %*% t(dXcentered), which is the covariance matrix with n x n.
whitenerXA = est.sigmaXA%^%(-0.5)   # ZCA Whitening, Lx. 
xDataA = whitenerXA %*% dXcentered   # Xw = Lx %*% Xc.matrix with n x px. 
invLx = est.sigmaXA%^%(0.5) # Inverse matrix of Lx, which is the whitenerXA aforemetioned. 

# For Y
# Scale rowwise
est.sigmaYA = tcrossprod(dYcentered)/(pY-1)  ## since already centered, can just take tcrossprod
whitenerYA = est.sigmaYA%^%(-0.5)   # ZCA Whitening
yDataA = whitenerYA %*% dYcentered   
invLy = est.sigmaYA%^%(0.5)


# Calculate JB values
JBall = calculateJB(matchMxMy$Ux[1:joint_rank, ], X = xDataA) + calculateJB(matchMxMy$Uy[1:joint_rank, ], X = yDataA) # the columns of Ux & Uy are up to the joint_rank

# SING on shared + individual
###################################
#R code for curvilinear search

# JB and tolerance parameters
alpha = 0.8
tol = 1e-10

# small rho
rho = JBall/10
out_indiv_small <- curvilinear(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = 2)
out_indiv_small_c <- SING::updateUboth_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = 2)
save(out_indiv_small, file = "out_indiv_small_example.Rda")


# small rho 
Sxmatched = out_indiv_small$Ux[1:2, ] %*% xDataA   # t() at here.
Symatched = out_indiv_small$Uy[1:2, ] %*% yDataA
Mxjoint = tcrossprod(invLx, matchMxMy$Ux[1:2, ])
Myjoint = tcrossprod(invLy, matchMxMy$Uy[1:2, ])


Sxtrue = data$sjX
Sytrue = data$sjY


### Signchange for the S matrix
Sxtrue = signchange(t(Sxtrue))
Sytrue = signchange(t(Sytrue))
Sx_rhoSmall = signchange(t(Sxmatched))
Sy_rhoSmall = signchange(t(Symatched))



lgrid = 33
par(mfrow = c(2,4))

image(matrix(Sxtrue[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth")
image(matrix(Sxtrue[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth")
image(vec2net(Sytrue[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sytrue[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth

image(matrix(Sx_rhoSmall[,1], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Small")
image(matrix(Sx_rhoSmall[,2], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="rho=Small")
image(vec2net(Sy_rhoSmall[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Small") 
image(vec2net(Sy_rhoSmall[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Small") 







