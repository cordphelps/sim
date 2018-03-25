simecol XY molecular diffusion (random walk)
================

``` r
library("simecol")
```

    ## Loading required package: deSolve

``` r
# https://cran.r-project.org/web/packages/simecol/simecol.pdf
diffusion.model <- rwalkModel(
  main = function(time, init, parms, inputs = NULL) {
    speed <- parms$speed
    xleft <- parms$area[1]
    xright <- parms$area[2]
    ybottom <- parms$area[3]
    ytop <- parms$area[4]
    x <- init$x # x coordinate
    y <- init$y # y coordinate
    a <- init$a # angle (in radians)
    n <- length(a)
    ## Rule 1: respect environment (grid as given in "inputs")
    ## 1a) identify location on "environmental 2D grid" for each individual
    i.j <- array(c(pmax(1, ceiling(x)), pmax(1, ceiling(y))), dim=c(n, 2))
    ## 1b) speed dependend on "environmental conditions"
    speed <- speed * inputs[i.j]
    ## Rule 2: Random Walk
    a <- (a + 2 * pi / runif(a))
    dx <- speed * cos(a)
    dy <- speed * sin(a)
    x <- x + dx
    y <- y + dy
    ## Rule 3: Wrap Around
    x <- ifelse(x > xright, xleft, x)
    y <- ifelse(y > ytop, ybottom, y)
    x <- ifelse(x < xleft, xright, x)
    y <- ifelse(y < ybottom, ytop, y)
    data.frame(x=x, y=y, a=a)
  },
  times = c(from=0, to=100, by=1),
  parms = list(ninds=50, speed = 1, area = c(0, 100, 0, 100)),
  solver = "iteration",
  initfunc = function(obj) {
    ninds <- obj@parms$ninds
    xleft <- obj@parms$area[1]
    xright <- obj@parms$area[2]
    ybottom <- obj@parms$area[3]
    ytop <- obj@parms$area[4]
    obj@init <- data.frame(x = runif(ninds) * (xright - xleft) + xleft,
                          y = runif(ninds) * (ytop - ybottom) + ybottom,
                          a = runif(ninds) * 2 * pi)
    inp <- matrix(1, nrow=100, ncol=100)
    inp[, 45:55] <- 0.2
    inputs(obj) <- inp
    obj
  }
)
```

``` r
diffusion.sim <- sim(diffusion.model)

out(diffusion.sim)[[34]]$x
```

    ##  [1] 89.01311736 17.35509104 38.89933667  8.61174157 38.56107661
    ##  [6] 24.71023626 30.02875948 48.39796311 65.34214774 48.18045579
    ## [11] 61.81210263 11.84914319 89.28195877 68.67182108 74.17385669
    ## [16]  6.17341074 21.32980205 56.98663007 61.20544184 63.13575690
    ## [21] 39.70584314 31.82293600 49.56178737 53.17552185 15.52038291
    ## [26] 80.12240186 66.83647257 94.01077226 58.04066604 20.13629380
    ## [31] 73.19891449  0.01089931 71.49542660 17.54140149 84.14100287
    ## [36] 80.67124307 98.82442989 17.79639744 71.48920893  8.03830982
    ## [41] 91.63653394 83.61048827  3.30958971 24.01592222 46.80807161
    ## [46] 28.71856715 43.73080681 94.42781953 59.26492166 69.62789049

``` r
plot(out(diffusion.sim)[[34]]$x)
```

![](mol-diffusion_files/figure-markdown_github/unnamed-chunk-3-1.png)

``` r
plotupca <- function(obj, type) {
  
  library(gridExtra)
  library(ggplot2)
  
  if (type=="S4") {
    o.df <- out(obj)  # output a dataframe
  } else {
    o.df <- obj
  }
  
  head(o.df)
  
  gg1 <- ggplot()  +
    geom_point(aes(y=rowNum, x=vectorX), o.df, colour = "purple", size=1, alpha=0.8)  +
    geom_point(aes(y=rowNum, x=vectorY), o.df, colour = "green", size=1, alpha=0.8)
  
  gg2 <- ggplot()  +
    geom_point(aes(x=vectorX, y=vectorY), o.df, fill="green", shape=21, colour = "purple", alpha=0.6) 
  
  
  grid.arrange(gg1, gg2, ncol=2, nrow=1)

}
```

``` r
# plotupca(out(diffusion.sim)[[34]]$x)

vectorX <- out(diffusion.sim)[[34]]$x
vectorY <- out(diffusion.sim)[[34]]$y

points.df <- as.data.frame(vectorX)
points.df$vectorY <- vectorY
points.df$rowNum <- seq.int(nrow(points.df))

head(points.df)
```

    ##     vectorX  vectorY rowNum
    ## 1 89.013117 23.25283      1
    ## 2 17.355091 70.10649      2
    ## 3 38.899337 30.02460      3
    ## 4  8.611742 44.15206      4
    ## 5 38.561077 57.58343      5
    ## 6 24.710236 63.53623      6

``` r
plotupca(points.df, type="S3")
```

![](mol-diffusion_files/figure-markdown_github/unnamed-chunk-5-1.png)

``` r
# plotupca(out(diffusion.sim)[[34]]$x)

# graph the path of a particle

# walk through the 50 moves (or "observations") 
# getting the x,y coords for each particle
move.list <- list()
for (j in 1:50) {
    vectorX <- out(diffusion.sim)[[j]]$x # get all x from observation j
    vectorY <- out(diffusion.sim)[[j]]$y
    path.df <- as.data.frame(vectorX)
    path.df$vectorY <- vectorY
    move.list[[j]] <- path.df
    rm("path.df")
}
# move.list is a list of dataframes that individually represent the x,y path
# of each particle (row1/df1 represents particle1)
# CAUTION: wrap-around is not identified

# plot the path of particle #21 and #41
df <- as.data.frame(move.list[21])
gg1 <- ggplot()  +
    geom_path(aes(x=vectorX, y=vectorY), df, colour = "purple", alpha=0.6) 
df <- as.data.frame(move.list[41])
gg2 <- ggplot()  +
    geom_path(aes(x=vectorX, y=vectorY), df, colour = "green", alpha=0.6)

grid.arrange(gg1, gg2, ncol=2, nrow=1)
```

![](mol-diffusion_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
head(out(diffusion.sim)[22])
```

    ## [[1]]
    ##             x         y         a
    ## 1  88.1069419 21.903629 1056.7034
    ## 2  18.8116680 70.491283  451.9840
    ## 3  37.0608377 27.347766  962.3657
    ## 4   5.7794797 42.783464  750.3946
    ## 5  42.7372491 59.137214  600.2553
    ## 6  25.7782510 65.773439  418.8945
    ## 7  30.3068019 44.110944 5544.4007
    ## 8  46.6616449 33.753496  673.5486
    ## 9  62.4325551 76.716339  614.8185
    ## 10 48.9068753  4.343983 2100.6617
    ## 11 61.8639800 59.043732  900.2822
    ## 12 13.5315581 62.032983  831.6806
    ## 13 88.6111528 36.904057  336.2852
    ## 14 68.5335120 93.779019  376.2736
    ## 15 72.3216858 80.039220 1812.7264
    ## 16  5.7362587 97.841703  501.1870
    ## 17 24.0556703 17.256601 2979.2857
    ## 18 56.0553789 88.174854  986.2690
    ## 19 62.8664621 19.309190  988.9566
    ## 20 63.8764422 92.931350  573.3307
    ## 21 36.3002415 81.811485 5548.2279
    ## 22 33.4412259 10.545159  452.7121
    ## 23 52.4866987 99.333618  423.6467
    ## 24 51.8468258 37.193486  389.6698
    ## 25 15.2769723 54.711893 7539.3866
    ## 26 82.1135812 14.946973  615.8095
    ## 27 63.4610771 75.121324 1065.8060
    ## 28 92.0779951 93.281529  615.2939
    ## 29 60.1973122 13.535993  389.8092
    ## 30 25.5008266 93.788673  468.5204
    ## 31 73.5055018 44.656011  534.5469
    ## 32  3.3261951 97.223264  251.8070
    ## 33 71.2306352 56.781904  534.2950
    ## 34 21.1428668 89.450619  475.8792
    ## 35 82.0159381 68.910794 6259.2100
    ## 36 81.5952713 17.260739 1828.1973
    ## 37  0.7626647 91.168229  446.9280
    ## 38 16.9320044 60.420545  471.3892
    ## 39 73.8808286 14.946369  496.5975
    ## 40  9.1189651  9.879275 1569.0033
    ## 41 94.7276022 99.922557  282.8154
    ## 42 82.9729114  8.082292  815.3122
    ## 43  4.1992440 38.882415  790.8475
    ## 44 24.9968887 32.972442 1962.5239
    ## 45 48.0448815 78.399923  850.2663
    ## 46 29.8272395 59.715871  634.3531
    ## 47 43.1425489 58.405969 7996.7735
    ## 48 94.3509172 52.758425 1060.3177
    ## 49 59.6085089 40.561467  833.8069
    ## 50 69.7803882 75.356244 1208.0552

``` r
main(diffusion.sim)
```

    ## function(time, init, parms, inputs = NULL) {
    ##     speed <- parms$speed
    ##     xleft <- parms$area[1]
    ##     xright <- parms$area[2]
    ##     ybottom <- parms$area[3]
    ##     ytop <- parms$area[4]
    ##     x <- init$x # x coordinate
    ##     y <- init$y # y coordinate
    ##     a <- init$a # angle (in radians)
    ##     n <- length(a)
    ##     ## Rule 1: respect environment (grid as given in "inputs")
    ##     ## 1a) identify location on "environmental 2D grid" for each individual
    ##     i.j <- array(c(pmax(1, ceiling(x)), pmax(1, ceiling(y))), dim=c(n, 2))
    ##     ## 1b) speed dependend on "environmental conditions"
    ##     speed <- speed * inputs[i.j]
    ##     ## Rule 2: Random Walk
    ##     a <- (a + 2 * pi / runif(a))
    ##     dx <- speed * cos(a)
    ##     dy <- speed * sin(a)
    ##     x <- x + dx
    ##     y <- y + dy
    ##     ## Rule 3: Wrap Around
    ##     x <- ifelse(x > xright, xleft, x)
    ##     y <- ifelse(y > ytop, ybottom, y)
    ##     x <- ifelse(x < xleft, xright, x)
    ##     y <- ifelse(y < ybottom, ytop, y)
    ##     data.frame(x=x, y=y, a=a)
    ##   }

``` r
init(diffusion.sim)
```

    ##             x          y           a
    ## 1  89.8176688 20.9175626 3.387127430
    ## 2  17.5608011 69.8094813 1.483246387
    ## 3  36.0628396 25.8613132 3.107248560
    ## 4   3.2248680 42.3110943 0.137712781
    ## 5  39.2654890 55.0047679 5.796740170
    ## 6  23.0106839 65.0140717 5.527337248
    ## 7  31.6099235 42.2819749 2.450001390
    ## 8  42.3162634 33.9713274 2.009350406
    ## 9  65.1827416 75.6699523 2.881626803
    ## 10 52.9802584  3.8122732 6.159252631
    ## 11 69.3161634 57.5289337 2.084791738
    ## 12 11.6383054 62.6657421 2.454157323
    ## 13 87.0053281 40.4591758 4.666168643
    ## 14 68.9865886 90.3023596 5.753449563
    ## 15 73.1338595 87.4366200 0.873787447
    ## 16  8.3815161 98.4746788 4.530808450
    ## 17 20.3488982 17.2032066 6.053773853
    ## 18 56.8749255 94.4561264 0.916619158
    ## 19 68.2660590 21.1796992 0.007671978
    ## 20 64.8897163 92.7621396 2.688439101
    ## 21 36.3935287 88.4673918 1.237992251
    ## 22 29.7668572  7.5557267 1.084638004
    ## 23 55.8973416  0.6793483 0.579682715
    ## 24 49.7417367 38.3464492 3.546079531
    ## 25 15.9644923 54.4672432 3.289748782
    ## 26 87.3064513 18.4615904 3.974448490
    ## 27 63.3028936 77.6207658 5.543305626
    ## 28 90.3343166 96.3508712 1.714354829
    ## 29 54.5171436 13.6325621 0.339602018
    ## 30 27.7560913  0.8946603 1.612901712
    ## 31 72.9767401 43.7358724 2.576238408
    ## 32  3.2180990 98.7733793 2.896746097
    ## 33 69.8418966 60.8247170 4.463884358
    ## 34 26.2880032 92.7181595 3.006795205
    ## 35 77.8352359 65.3100064 0.535510467
    ## 36 78.8171917 20.0019937 3.225144307
    ## 37 98.2964007 91.2003888 4.075325148
    ## 38 10.5916161 56.2403478 3.098800737
    ## 39 70.2099336 18.9800677 4.978664838
    ## 40  8.7961914  8.0501889 5.720089552
    ## 41  0.8163733  0.9506251 5.425969913
    ## 42 89.5824264 11.6603815 2.413162890
    ## 43  2.3256850 39.3538766 1.998334728
    ## 44 22.7789989 36.0522785 4.067705643
    ## 45 50.3753918 73.0035759 2.695442738
    ## 46 28.1467442 63.3978484 2.195689804
    ## 47 44.2778106 63.7572796 2.474120815
    ## 48 93.9510361 52.3763277 0.274196285
    ## 49 62.1560400 42.5575109 5.380322927
    ## 50 69.9481701 73.9851655 2.780082877

``` r
parms(diffusion.sim)
```

    ## $ninds
    ## [1] 50
    ## 
    ## $speed
    ## [1] 1
    ## 
    ## $area
    ## [1]   0 100   0 100

``` r
equations(diffusion.sim)
```

    ## NULL

``` r
solver(diffusion.sim)
```

    ## [1] "iteration"

``` r
class(diffusion.sim)
```

    ## [1] "rwalkModel"
    ## attr(,"package")
    ## [1] "simecol"

``` r
#str(diffusion.sim)
```
