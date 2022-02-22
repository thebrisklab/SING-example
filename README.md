SING example tutorial
================
Liangkang Wang
12/06/2021

## Introduction

SING is used to extract joint and individual non-gaussian components
from different datasets. This is a tutorial example supporting the paper
**Simultaneous Non-Gaussian Component Analysis (SING) for Data
Integration in Neuroimaging Benjamin Risk, Irina Gaynanova**
<https://arxiv.org/abs/2005.00597v1>

The function codes are saved in `jngcaFunctions.R` and
`generate_data.R`. The data processing is saved in `My example.R` while
the figures of outcome are saved in `My example2.R`.

## Example Code

### Generate data

``` r
source("jngcaFunctions.R")
source('generate_data.R')
set.seed(0573452)
data <- generateData_v3(nsubject = 48, snr = c(1, 1), vars = c(0.005, 0.005))

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
```

![](figs/Original%20data%20figure.png)<!-- -->

### Pipeline of SING method

The data process code can be found in the two code file mentioned
before.

### Outcome

#### Joint X components

``` r
source("jngcaFunctions.R")
load("EstimatedComponents_example.Rda")
load("Estimated joint subject score_example.Rda")


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
```

![](figs/Joint%20X%20component.png)<!-- -->

### Joint Y components

``` r
par(mfrow = c(2,4))
image(vec2net(Sytrue[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sytrue[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="Truth") #Truth
image(vec2net(Sy_rhoSmall[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Large") # large rho
image(vec2net(Sy_rhoSmall[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="ro0=Large") # large rho
image(vec2net(SyjointICA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA") # Joint ICA
image(vec2net(SyjointICA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="jointICA") #Joint ICA
image(vec2net(SymCCA[,1]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA") #mCCA+jICA
image(vec2net(SymCCA[,2]), col = heat.colors(12), xaxt = "n", yaxt = "n",main="mCCA") #mCCA+jICA
```

![](figs/Joint%20Y%20component.png)<!-- -->

### Subject Scores

``` r
###Figure for the joint subject score


library(tidyverse)
library(ggpubr)

#True mj
t1 <- ggplot(data = trueMj)+
  geom_point(mapping = aes(y=mj1,x=number))+
  ggtitle("TrueMj,Comp1")+
  theme_bw()+
  theme(panel.grid = element_blank())

t2 <- ggplot(data = trueMj)+
  geom_point(mapping = aes(y=mj2,x=number))+
  ggtitle("TrueMj,Comp2")+
  theme_bw()+
  theme(panel.grid = element_blank())

#SING mj
S1 <- ggplot(data = SINGMj)+
  geom_point(mapping = aes(y=mj1,x=number))+
  ggtitle("SINGMj,Comp1")+
  theme_bw()+
  theme(panel.grid = element_blank())

S2 <- ggplot(data = SINGMj)+
  geom_point(mapping = aes(y=mj2,x=number))+
  ggtitle("SINGmj,Comp2")+
  theme_bw()+
  theme(panel.grid = element_blank())

#ICA mj
I1 <- ggplot(data = ICAMj)+
  geom_point(mapping = aes(y=mj2,x=number))+
  ggtitle("Joint Score,Comp1")+
  theme_bw()+
  theme(panel.grid = element_blank())  

#due to the permutation test, the sequence at here is inverse, which can be seen in the previous figure of joint components

I2 <- ggplot(data = ICAMj)+
  geom_point(mapping = aes(y=mj1,x=number))+
  ggtitle("Joint Score,Comp2")+
  theme_bw()+
  theme(panel.grid = element_blank())

ggarrange(t1,S1,I1,t2,S2,I2,ncol = 3,nrow = 2)
```

![](figs/Subject%20Score.png)<!-- -->
