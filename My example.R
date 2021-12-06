# Generate simulation data

# Source all the functions
source("jngcaFunctions.R")
source("mCCAjointICA.R")
source('generate_data.R')
# Generate data, will return both components and mixing matrix
set.seed(0573452)
data <- generateData_v3(nsubject = 48, snr = c(1, 1), vars = c(0.005, 0.005))
save(data,file='data_example.Rda')
# Components for X
###simData

# simData=newSimFMRI()

#simData=newSimFMRI()


#par(mfrow = c(1,4))
#image(matrix(simData$S[,1],33))
#image(matrix(simData$S[,2],33))
#image(matrix(simData$S[,3],33))
#image(matrix(simData$S[,4],33))

lgrid = 33
par(mfrow = c(2,4))
# Components for X
image(matrix(data$sjX[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="X joint component1")
image(matrix(data$sjX[2,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="X joint component2")
image(matrix(data$siX[1,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="X independent component2")
image(matrix(data$siX[2,], lgrid, lgrid), col = heat.colors(12), xaxt = "n", yaxt = "n",main="X independent component2")

# Components for Y
image(vec2net(data$sjY[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Y joint component1")
image(vec2net(data$sjY[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Y joint component2")
image(vec2net(data$siY[1,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Y independent component1")
image(vec2net(data$siY[2,]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Y independent component2")


# Use SING to conmpute
# Center X and Y
n = nrow(data$dX)
pX = ncol(data$dX)
pY = ncol(data$dY)
dXcentered <- data$dX - matrix(rowMeans(data$dX), n, pX, byrow = F)
dYcentered <- data$dY - matrix(rowMeans(data$dY), n, pY, byrow = F)
trueJx <- data$mj %*% data$sjX
trueJy <- t(t(data$mj) * c(-5, 2)) %*% data$sjY
trueJxF2 <- sum(trueJx^2)
trueJyF2 <- sum(trueJy^2)

## Apply JointICA
#############################################################
out_jointICA <- jointICA(dXcentered, dYcentered, r0 = 2)

# Prepare the output list for joint ICA
save(out_jointICA, file = "jointICA_example.Rda")

# Joint rank estimation based on oversaturated models
#####################################################

estX_JB = mlcaFP(xData = t(data$dX), n.comp = 48, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')

Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(data$dX))

estY_JB = mlcaFP(xData = t(data$dY), n.comp = 48, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')

My_JB = est.M.ols(sData = estY_JB$S, xData = t(data$dY))

# Get joint components out
alpha = 0.05
nperms = 10 #1000
matchedResults = angleMatchICA(t(Mx_JB), t(My_JB))
permuteJoint = permmatRank_joint(matchedResults, nperms = nperms)
joint_rank = min(which(permuteJoint$pvalues > alpha)) - 1
pval_joint = permuteJoint$pvalues
joint_rank # selects rank 2

## Apply separate JB 
#############################################################
# JB on X
estX_JB = mlcaFP(xData = t(data$dX), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
Uxfull <- estX_JB$Ws

# Match mixing matrices to get correct ordering, then can get starting points
Mx_JB = est.M.ols(sData = estX_JB$S, xData = t(data$dX))

# JB on Y
estY_JB = mlcaFP(xData = t(data$dY), n.comp = 12, whiten = 'sqrtprec', restarts.pbyd = 20, distribution='JB')
Uyfull <- estY_JB$Ws

# Get joint components out
My_JB = est.M.ols(sData = estY_JB$S, xData = t(data$dY))

matchMxMy = greedymatch(t(Mx_JB), t(My_JB), Ux = t(Uxfull), Uy = t(Uyfull))

cor(matchMxMy$Mx[, 1:2], matchMxMy$My[, 1:2]) # 0.95 and 0.93

errorMx = frobICA(M1 = t(matchMxMy$Mx[, 1:2]), M2 = t(data$mj), standardize = T)*sqrt(ncol(Mx_JB))
errorMx # 0.295

errorMy = frobICA(M1 = t(matchMxMy$My[, 1:2]), M2 = t(data$mj), standardize = T)*sqrt(ncol(My_JB))
errorMy # 0.208

errorSx = frobICA(S1 = t(data$sjX), S2 = estX_JB$S[, matchMxMy$mapX[1:2]], standardize  = T)
errorSx # 0.051 

errorSy = frobICA(S1 = t(data$sjY), S2 = estY_JB$S[, matchMxMy$mapY[1:2]], standardize  = T)
errorSy # 0.031

errorJx = sum((tcrossprod(matchMxMy$Mx[, 1:2],  estX_JB$S[, matchMxMy$mapX[1:2]]) - trueJx)^2)/trueJxF2
errorJx #0.095

errorJy = sum((tcrossprod(matchMxMy$My[, 1:2], estY_JB$S[, matchMxMy$mapY[1:2]])  - trueJy)^2)/trueJyF2
errorJy #0.025


# Prepare the output list for separate JB
output_sepJB <- list(errorSx = errorSx, errorSy = errorSy, errorMx = errorMx, errorMy = errorMy, errorJx = errorJx, errorJy = errorJy, estX_JB = estX_JB, estY_JB = estY_JB, mj = data$mj)
save(output_sepJB, file = "sepJB_example.Rda")


# Joint rank estimation based on separate models fitted with true number of components
#####################################################
# Rank estimation via permutation
alpha = 0.05
nperms = 10
matchedResults = angleMatchICA(t(Mx_JB), t(My_JB))
permuteJoint = permmatRank_joint(matchedResults, nperms = nperms)
joint_rank = min(which(permuteJoint$pvalues > alpha)) - 1
pval_joint = permuteJoint$pvalues
joint_rank # selects rank 2


# Whiten dX and dY
###################################

# Scale rowwise
est.sigmaXA = tcrossprod(dXcentered)/(pX-1)  ## since already centered, can just take tcrossprod
whitenerXA = est.sigmaXA%^%(-0.5)
invLx = est.sigmaXA%^%(0.5)
xDataA = whitenerXA %*% dXcentered

# For Y
# Scale rowwise
est.sigmaYA = tcrossprod(dYcentered)/(pY-1)  ## since already centered, can just take tcrossprod
whitenerYA = est.sigmaYA%^%(-0.5)
yDataA = whitenerYA %*% dYcentered
invLy = est.sigmaYA%^%(0.5)

# Starting points for the algorithm are Us
##########################################

# Calculate JB values
JBall = calculateJB(matchMxMy$Ux[1:2, ], X = xDataA) + calculateJB(matchMxMy$Uy[1:2, ], X = yDataA)


# SING on shared + individual
###################################
# Libraries needed to run C functions
library(Rcpp)
library(RcppArmadillo)

sourceCpp("CfunctionsCurvilinear.cpp")

# JB and tolerance parameters
alpha = 0.8
tol = 1e-10

# small rho
rho = JBall/10
out_indiv_small <- updateUboth_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = matchMxMy$Ux, Uy = matchMxMy$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = 2)
save(out_indiv_small, file = "out_indiv_small_example.Rda")

# medium rho
rho = JBall
out_indiv_medium <- updateUboth_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = out_indiv_small$Ux, Uy = out_indiv_small$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = 2)
save(out_indiv_medium, file = "out_indiv_medium_example.Rda")

# large rho
rho = JBall * 10
out_indiv_large <- updateUboth_c(invLx = invLx, invLy = invLy, xData = xDataA, yData = yDataA, Ux = out_indiv_medium$Ux, Uy = out_indiv_medium$Uy, rho = rho, tol = tol, alpha = alpha, maxiter = 1500, r0 = 2)
save(out_indiv_large, file = "out_indiv_large_example.Rda")



###############################################################################
# Apply mCCA + Joint ICA, and calculate the errors
###############################################################################
out_mcca <- mCCA_jointICA(dXcentered, dYcentered, Mx = 12, My = 12, M = 2)

save(out_mcca, file = "mCCA_LargeScale_example.Rda")



# Error calculation for all methods, and creation of figures

# Load the results of joint ICA
load("jointICA_example.Rda")
# Load the results from separate approach (rho = 0)
load("sepJB_example.Rda")
# Load the results from SING (from small rho to large rho)
load("out_indiv_small_example.Rda")
load("out_indiv_medium_example.Rda")
load("out_indiv_large_example.Rda")
# Load the results from mCCA + joint ICA
load("mCCA_LargeScale_example.Rda")

# Load true underlying components
load("data_example.Rda")
trueJx <- data$mj %*% data$sjX
trueJy <- t(t(data$mj) * c(-5, 2)) %*% data$sjY
trueJxF2 <- sum(trueJx^2)
trueJyF2 <- sum(trueJy^2)

# Source all the functions
source("jngcaFunctions.R")
source("mCCAjointICA.R")


# Create data frame to store results for all methods
methods = c("Joint ICA", "mCCA+jICA", "rho 0", "rho averaged", "small rho", "medum rho", "large rho")
types = c("Sx", "Sy", "Mx", "My", "Jx", "Jy")
errors = matrix(NA, length(methods), length(types))
errors = as.data.frame(errors)
colnames(errors) = types
rownames(errors) = methods

# Data processing and dimensions
n = nrow(data$dX)
pX = ncol(data$dX)
pY = ncol(data$dY)
# For X
# Center
dXsA <- data$dX - matrix(rowMeans(data$dX), n, pX, byrow = F)
# Scale rowwise
est.sigmaXA = tcrossprod(dXsA)/(pX-1)  ## since already centered, can just take tcrossprod
whitenerXA = est.sigmaXA%^%(-0.5)
invLx = est.sigmaXA%^%(0.5)
xDataA = whitenerXA %*% dXsA

# For Y
# Center
dYsA <- data$dY - matrix(rowMeans(data$dY), n, pY, byrow = F)
# Scale rowwise
est.sigmaYA = tcrossprod(dYsA)/(pY-1)  ## since already centered, can just take tcrossprod
whitenerYA = est.sigmaYA%^%(-0.5)
yDataA = whitenerYA %*% dYsA
invLy = est.sigmaYA%^%(0.5)


# Starting points for the algorithm are Us
###########################################
Uxfull <- output_sepJB$estX_JB$Ws
Mx_JB = est.M.ols(sData = output_sepJB$estX_JB$S, xData = t(data$dX))
Uyfull <- output_sepJB$estY_JB$Ws
My_JB = est.M.ols(sData = output_sepJB$estY_JB$S, xData = t(data$dY))

matchMxMy = greedymatch(t(Mx_JB), t(My_JB), Ux = t(Uxfull), Uy = t(Uyfull))

cor(matchMxMy$Mx[, 1:2], matchMxMy$My[, 1:2]) # 0.95 and 0.93

# Calculate errors for Sx
############################
#joint ICA
errorSx = frobICA(S1 = t(data$sjX), S2 = out_jointICA$S[1:pX, ], standardize  = T)
errorSx # 0.375
errors[1, 1] = errorSx

# mCCA + joint ICA
errorSx = frobICA(S1 = t(data$sjX), S2 = out_mcca$S[1:pX, ], standardize  = T)
errorSx #0.38
errors[2, 1] = errorSx

# rho = 0 with all 2 matched components
Sxmatched = matchMxMy$Ux[1:2, ] %*% xDataA
errorSx = frobICA(S1 = t(data$sjX), S2 = t(Sxmatched), standardize  = T)
errorSx # 0.178
errors[3, 1] = errorSx

# small rho
Sxmatched = out_indiv_small$Ux[1:2, ] %*% xDataA
errorSx = frobICA(S1 = t(data$sjX), S2 = t(Sxmatched), standardize  = T)
errorSx # 0.176
errors[5, 1] = errorSx

# medium rho
Sxmatched = out_indiv_medium$Ux[1:2, ] %*% xDataA
errorSx = frobICA(S1 = t(data$sjX), S2 = t(Sxmatched), standardize  = T)
errorSx # 0.163
errors[6, 1] = errorSx

# large rho
Sxmatched = out_indiv_large$Ux[1:2, ] %*% xDataA
errorSx = frobICA(S1 = t(data$sjX), S2 = t(Sxmatched), standardize  = T)
errorSx # 0.16
errors[7, 1] = errorSx

# Calculate errors for Sy
############################
#joint ICA
errorSy = frobICA(S1 = t(data$sjY), S2 = out_jointICA$S[(pX+1):(pX+pY), ], standardize  = T)
errorSy # 0.399
errors[1, 2] = errorSy

# mCCA + joint ICA
errorSy = frobICA(S1 = t(data$sjY), S2 = out_mcca$S[(pX+1):(pX+pY), ], standardize  = T)
errorSy #0.402
errors[2, 2] = errorSy

# rho = 0 with all 2 matched components
Symatched = matchMxMy$Uy[1:2, ] %*% yDataA
errorSy = frobICA(S1 = t(data$sjY), S2 = t(Symatched), standardize  = T)
errorSy # 0.059
errors[3, 2] = errorSy

# small rho 
Symatched = out_indiv_small$Uy[1:2, ] %*% yDataA
errorSy = frobICA(S1 = t(data$sjY), S2 = t(Symatched), standardize  = T)
errorSy # 0.059
errors[5, 2] = errorSy


# medium rho
Symatched = out_indiv_medium$Uy[1:2, ] %*% yDataA
errorSy = frobICA(S1 = t(data$sjY), S2 = t(Symatched), standardize  = T)
errorSy # 0.067
errors[6, 2] = errorSy

# large rho
Symatched = out_indiv_large$Uy[1:2, ] %*% yDataA
errorSy = frobICA(S1 = t(data$sjY), S2 = t(Symatched), standardize  = T)
errorSy # 0.085
errors[7, 2] = errorSy


# Calculate errors for M
############################
# joint ICA
dXcentered <- data$dX - matrix(rowMeans(data$dX), n, pX, byrow = F)
dYcentered <- data$dY - matrix(rowMeans(data$dY), n, pY, byrow = F)
# Normalization step here. Divide by \|X\|_F^2 each part
dXsA <- dXcentered/sqrt(mean(dXcentered^2))
dYsA <- dYcentered/sqrt(mean(dYcentered^2))
# Concatenate them together [X, Y] and perform SVD (PCA)
dXY <- cbind(dXsA, dYsA) # [X, Y] ~ UDV'
Mjoint = est.M.ols(sData = out_jointICA$S, xData = t(dXY))
errorM = frobICA(M1 = Mjoint, M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.123
errors[1, 3:4] = errorM

# mCCA + joint ICA
errorMx = frobICA(M1 = out_mcca$Mx, M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
errorMx #0.178
errors[2, 3] = errorMx
errorMy = frobICA(M1 = out_mcca$My, M2 = t(data$mj), standardize  = T) * sqrt(nrow(data$mj))
errorMy #0.156
errors[2, 4] = errorMy

# rho = 0 with on matched
Mxjoint = tcrossprod(invLx, matchMxMy$Ux[1:2, ])
errorMx = frobICA(M1 = t(Mxjoint), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorMx # 0.0643

Myjoint = tcrossprod(invLy, matchMxMy$Uy[1:2, ])
errorMy = frobICA(M1 = t(Myjoint), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorMy # 0.034

errors[3, 3] = errorMx
errors[3, 4] = errorMy

cordiag = diag(cor(Mxjoint, Myjoint))
cordiag # 0.95 and 0.93
newM = (Mxjoint + Myjoint%*% sign(diag(cordiag)))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.029

# small rho on matched
Mxjoint = tcrossprod(invLx, out_indiv_small$Ux[1:2, ])
errorMx = frobICA(M1 = t(Mxjoint), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorMx # 0.0614

Myjoint = tcrossprod(invLy, out_indiv_small$Uy[1:2, ])
errorMy = frobICA(M1 = t(Myjoint), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorMy # 0.031

errors[5, 3] = errorMx
errors[5, 4] = errorMy

cordiag = diag(cor(Mxjoint, Myjoint))
cordiag #0.9986; 0.9969
newM = (Mxjoint + Myjoint%*% sign(diag(cordiag)))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.028

# medium rho on matched
Mxjoint = tcrossprod(invLx, out_indiv_medium$Ux[1:2, ])
errorMx = frobICA(M1 = t(Mxjoint), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorMx # 0.0482

Myjoint = tcrossprod(invLy, out_indiv_medium$Uy[1:2, ])
errorMy = frobICA(M1 = t(Myjoint), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorMy # 0.024

errors[6, 3] = errorMx
errors[6, 4] = errorMy

cordiag = diag(cor(Mxjoint, Myjoint))
cordiag #0.995895; 0.998219
newM = (Mxjoint + Myjoint%*% sign(diag(cordiag)))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.0243

# large rho on matched
Mxjoint = tcrossprod(invLx, out_indiv_large$Ux[1:2, ])
errorMx = frobICA(M1 = t(Mxjoint), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorMx # 0.0278

Myjoint = tcrossprod(invLy, out_indiv_large$Uy[1:2, ])
errorMy = frobICA(M1 = t(Myjoint), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorMy # 0.0225

errors[7, 3] = errorMx
errors[7, 4] = errorMy

cordiag = diag(cor(Mxjoint, Myjoint))
cordiag #0.999999999; 1
newM = (Mxjoint + Myjoint%*% sign(diag(cordiag)))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
newM = newM %*% diag(1/apply(newM, 2, function(x) sqrt(sum(x^2))))
errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize = T) * sqrt(ncol(Mjoint))
errorM # 0.0229



# Calculate errors for Jx and Jy
############################
# Joint ICA
errorJx = sum((tcrossprod(t(out_jointICA$Mjoint), out_jointICA$S[1:pX, ]) * sqrt(mean(dXcentered^2)) - trueJx)^2)/trueJxF2
errorJx #0.121
errors[1, 5] = sqrt(errorJx)
errorJy = sum((tcrossprod(t(out_jointICA$Mjoint), out_jointICA$S[(pX + 1):(pX + pY), ])* sqrt(mean(dYcentered^2)) - trueJy)^2)/trueJyF2 
errorJy # 0.081
errors[1, 6] = sqrt(errorJy)

# mCCA + Joint ICA
errorJx = sum((tcrossprod(t(out_mcca$Mx), out_mcca$S[1:pX, ]) - trueJx)^2)/trueJxF2
errorJx # 0.145
errors[2, 5] = sqrt(errorJx)
errorJy = sum((tcrossprod(t(out_mcca$My), out_mcca$S[(pX+1):(pX+pY), ]) - trueJy)^2)/trueJyF2
errorJy # 0.083 
errors[2, 6] = sqrt(errorJy)


# separate rho
Sxmatched = matchMxMy$Ux[1:2, ] %*% xDataA
Symatched = matchMxMy$Uy[1:2, ] %*% yDataA
Mxjoint = tcrossprod(invLx, matchMxMy$Ux[1:2, ])
Myjoint = tcrossprod(invLy, matchMxMy$Uy[1:2, ])

errorJx = sum((Mxjoint %*% Sxmatched  - trueJx)^2)/trueJxF2
errorJx #0.034
errorJy = sum((Myjoint %*% Symatched - trueJy)^2)/trueJyF2 
errorJy # 0.00554
errors[3, 5] = sqrt(errorJx)
errors[3, 6] = sqrt(errorJy)


# averaged rho
newM = aveM(matchMxMy$Mx[,1:2], matchMxMy$My[,1:2])
errorM = frobICA(M1 = t(newM), M2 = t(data$mj), standardize  = T)*sqrt(ncol(Mjoint))

errors[4, 3:4] = errorM

# Back projection on the mixing matrix
outx <- est.S.backproject(Sxmatched, Mxjoint, newM)
outy <- est.S.backproject(Symatched, Myjoint, newM)

errorSxAve = frobICA(S2 = t(data$sjX), S1 = t(outx$S), standardize  = T)
errorSyAve = frobICA(S2 = t(data$sjY), S1 = t(outy$S), standardize  = T)

errors[4, 1] = errorSxAve
errors[4, 2] = errorSyAve

# Joint signal Frobenius norm reconstruction error with average
errorJxave = sum((newM %*% diag(outx$D) %*% outx$S/sqrt(pX-1) - trueJx)^2)/trueJxF2
errorJyave = sum((newM %*% diag(outy$D) %*% outy$S/sqrt(pY-1) - trueJy)^2)/trueJyF2


errors[4, 5] = sqrt(errorJxave)
errors[4, 6] = sqrt(errorJyave)

# small rho
Sxmatched = out_indiv_small$Ux[1:2, ] %*% xDataA
Symatched = out_indiv_small$Uy[1:2, ] %*% yDataA
Mxjoint = tcrossprod(invLx, out_indiv_small$Ux[1:2, ])
Myjoint = tcrossprod(invLy, out_indiv_small$Uy[1:2, ])

errorJx = sum((Mxjoint %*% Sxmatched  - trueJx)^2)/trueJxF2
errorJx #0.0334
errorJy = sum((Myjoint %*% Symatched - trueJy)^2)/trueJyF2 
errorJy # 0.0055
errors[5, 5] = sqrt(errorJx)
errors[5, 6] = sqrt(errorJy)

# medium rho
Sxmatched = out_indiv_medium$Ux[1:2, ] %*% xDataA
Symatched = out_indiv_medium$Uy[1:2, ] %*% yDataA
Mxjoint = tcrossprod(invLx, out_indiv_medium$Ux[1:2, ])
Myjoint = tcrossprod(invLy, out_indiv_medium$Uy[1:2, ])

errorJx = sum((Mxjoint %*% Sxmatched  - trueJx)^2)/trueJxF2
errorJx #0.0283
errorJy = sum((Myjoint %*% Symatched - trueJy)^2)/trueJyF2 
errorJy # 0.0056
errors[6, 5] = sqrt(errorJx)
errors[6, 6] = sqrt(errorJy)

# large rho
Sxmatched = out_indiv_large$Ux[1:2, ] %*% xDataA
Symatched = out_indiv_large$Uy[1:2, ] %*% yDataA
Mxjoint = tcrossprod(invLx, out_indiv_large$Ux[1:2, ])
Myjoint = tcrossprod(invLy, out_indiv_large$Uy[1:2, ])

errorJx = sum((Mxjoint %*% Sxmatched  - trueJx)^2)/trueJxF2
errorJx #0.027
errorJy = sum((Myjoint %*% Symatched - trueJy)^2)/trueJyF2 
errorJy # 0.0064
errors[7, 5] = sqrt(errorJx)
errors[7, 6] = sqrt(errorJy)

errors = round(errors, 3)

# Create a table with all errors
####################################
library(xtable)
table = xtable(errors, digits = 3)
print(table)

# Save true components, estimated from jointICA, estimated with rho=0, estimated with small rho'
############################
Sxtrue = t(data$sjX)
Sytrue = t(data$sjY)
SxjointICA = out_jointICA$S[1:pX, ]
SyjointICA = out_jointICA$S[(pX+1):(pX + pY), ]
Sx_rho0 = t(matchMxMy$Ux[1:2, ] %*% xDataA)
Sy_rho0 = t(matchMxMy$Uy[1:2, ] %*% yDataA)
Sx_rhoSmall = t(out_indiv_small$Ux[1:2, ] %*% xDataA)  # I have make change at here to change the sequence of graph
Sy_rhoSmall = t(out_indiv_small$Uy[1:2, ] %*% yDataA)
Sx_rhoLarge = t(out_indiv_large$Ux[1:2, ] %*% xDataA)
Sy_rhoLarge = t(out_indiv_large$Uy[1:2, ] %*% yDataA)
SxmCCA = out_mcca$S[1:pX, ]
SymCCA = out_mcca$S[(pX+1):(pX + pY), ]

Sxtrue = signchange(Sxtrue)
Sytrue = signchange(Sytrue)
SxjointICA = signchange(SxjointICA)
SyjointICA = signchange(SyjointICA)
Sx_rhoSmall = signchange(Sx_rhoSmall)
Sy_rhoSmall = signchange(Sy_rhoSmall)
Sx_rhoLarge  = signchange(Sx_rhoLarge)
Sy_rhoLarge  = signchange(Sy_rhoLarge)
Sx_rho0  = signchange(Sx_rho0)
Sy_rho0  = signchange(Sy_rho0)
SxmCCA = signchange(SxmCCA)
SymCCA = signchange(SymCCA)

# Save the components
save(Sxtrue, Sytrue, SxjointICA, SyjointICA, Sx_rho0, Sy_rho0, Sx_rhoSmall, Sy_rhoSmall, Sx_rhoLarge, Sy_rhoLarge, SxmCCA, SymCCA,file = "EstimatedComponents_example.Rda")


## joint subject score
trueMj <- data.frame(mj1=data$mj[,1],mj2=data$mj[,2],number=1:48)
ICAMj <- data.frame(mj1=out_jointICA$Mjoint[1,],mj2=out_jointICA$Mjoint[2,],number=1:48)
SINGMj <- data.frame(mj1=newM[,1],mj2=newM[,2],number=1:48) # averaged rho

save(trueMj=trueMj,ICAMj=ICAMj,SINGMj=SINGMj,file = "Estimated joint subject score_example.Rda")






