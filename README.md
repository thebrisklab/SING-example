SING example tutorial
================
Liangkang Wang
12/06/2021

## Introduction

SING method is used in extracting joint non-gaussian components from
different datasets when wiping off gaussian noise. This is a tutorial
example supporting the paper **Simultaneous Non-Gaussian Component
Analysis (SING) for Data Integration in Neuroimaging Benjamin Risk,
Irina Gaynanova** <https://arxiv.org/abs/2005.00597v1>

The function codes are saved in `jngcaFunctions.R` and
`generate_data.R`. The data processing is saved in `My example.R` while
the figures of outcome are saved in My example2.R\`.

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

### Data process

The data process code can be found in the two code file mentioned
before.
