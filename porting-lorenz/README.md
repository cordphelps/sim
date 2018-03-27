simecol lorenz
================

### vensim mode 
 

![alt text](lorenz/Lorenz.png "optional title")


### simecol mode 

``` r
library("simecol")
```

    ## Loading required package: deSolve

``` r
lorenz.model <- new("odeModel",                   # object class = 'odeModel'
                
         main = function(times, y, parms) {   # object method
                                              # simecol magic: 2nd argument 'y' points to 'init'
          with(as.list(c(parms, y)), {

            horizTempDt <- ( convectiveFlowDelta * ( rayleigh - verticalTempDelta ) ) - horizTempDelta
            verticalTempDt <- (convectiveFlowDelta * horizTempDelta) - (height * verticalTempDelta)
            convectiveFlowDt <- prandtl * ( horizTempDelta - convectiveFlowDelta )
            
            horizTempDelta <- horizTempDt
            verticalTempDelta <- verticalTempDt
            convectiveFlowDelta <- convectiveFlowDt
          
            list(c(horizTempDelta, verticalTempDelta, convectiveFlowDelta))
           })
         },
         
         #.0078125
        times = seq(0, 50, .0078125), 
        parms = c(prandtl=10, rayleigh=28, height=8/3),
        init = c(convectiveFlowDelta=0, horizTempDelta=1, verticalTempDelta=0),
        solver = "rk4"
        
      )
```

``` r
lorenz.sim <- sim(lorenz.model)
#plot(lorenz.sim)
```

``` r
plotupca <- function(obj, ...) {
  
  library(gridExtra)
  library(ggplot2)
  
  o.df <- out(obj)  # output a dataframe
  
  gg1 <- ggplot()  +
    geom_line(aes(x=verticalTempDelta, y=convectiveFlowDelta), o.df, colour = "green", size=1, alpha=0.8)
  
  gg2 <- ggplot()  +
    geom_line(aes(x=time, y=convectiveFlowDelta), o.df,  colour = "purple", alpha=0.6) 
  
  
  grid.arrange(gg1, gg2, ncol=2, nrow=1)

}
```

``` r
 plotupca(lorenz.sim)
```

    ## Warning: Removed 6301 rows containing missing values (geom_path).

    ## Warning: Removed 6301 rows containing missing values (geom_path).

![](lorenz_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
head(out(lorenz.sim))
```

    ##        time convectiveFlowDelta horizTempDelta verticalTempDelta
    ## 1 0.0000000         0.000000000      1.0000000         0.0000000
    ## 2 0.0078125        -0.008728548      0.9991512         0.0784315
    ## 3 0.0156250        -0.019563975      0.9965861         0.1574912
    ## 4 0.0234375        -0.032988830      0.9922734         0.2372266
    ## 5 0.0312500        -0.049596997      0.9861773         0.3177246
    ## 6 0.0390625        -0.070118823      0.9782553         0.3991204

``` r
main(lorenz.sim)
```

    ## function(times, y, parms) {   # object method
    ##                                               # simecol magic: 2nd argument 'y' points to 'init'
    ##           with(as.list(c(parms, y)), {
    ## 
    ##             horizTempDt <- ( convectiveFlowDelta * ( rayleigh - verticalTempDelta ) ) - horizTempDelta
    ##             verticalTempDt <- (convectiveFlowDelta * horizTempDelta) - (height * verticalTempDelta)
    ##             convectiveFlowDt <- prandtl * ( horizTempDelta - convectiveFlowDelta )
    ##             
    ##             horizTempDelta <- horizTempDt
    ##             verticalTempDelta <- verticalTempDt
    ##             convectiveFlowDelta <- convectiveFlowDt
    ##           
    ##             list(c(horizTempDelta, verticalTempDelta, convectiveFlowDelta))
    ##            })
    ##          }

``` r
init(lorenz.sim)
```

    ## convectiveFlowDelta      horizTempDelta   verticalTempDelta 
    ##                   0                   1                   0

``` r
parms(lorenz.sim)
```

    ##   prandtl  rayleigh    height 
    ## 10.000000 28.000000  2.666667

``` r
equations(lorenz.sim)
```

    ## NULL

``` r
solver(lorenz.sim)
```

    ## [1] "rk4"

``` r
class(lorenz.sim)
```

    ## [1] "odeModel"
    ## attr(,"package")
    ## [1] "simecol"

``` r
str(lorenz.sim)
```

    ## Formal class 'odeModel' [package "simecol"] with 10 slots
    ##   ..@ parms    : Named num [1:3] 10 28 2.67
    ##   .. ..- attr(*, "names")= chr [1:3] "prandtl" "rayleigh" "height"
    ##   ..@ init     : Named num [1:3] 0 1 0
    ##   .. ..- attr(*, "names")= chr [1:3] "convectiveFlowDelta" "horizTempDelta" "verticalTempDelta"
    ##   ..@ observer : NULL
    ##   ..@ main     :function (times, y, parms)  
    ##   .. ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 6 17 20 10 17 10 6 20
    ##   .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x7f8bf2d29e98> 
    ##   ..@ equations: NULL
    ##   ..@ times    : num [1:6401] 0 0.00781 0.01562 0.02344 0.03125 ...
    ##   ..@ inputs   : NULL
    ##   ..@ solver   : chr "rk4"
    ##   ..@ out      : deSolve [1:6401, 1:4] 0 0.00781 0.01562 0.02344 0.03125 ...
    ##   .. ..- attr(*, "istate")= int [1:21] 0 6400 25601 NA NA NA NA NA NA NA ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   .. .. ..$ : NULL
    ##   .. .. ..$ : chr [1:4] "time" "convectiveFlowDelta" "horizTempDelta" "verticalTempDelta"
    ##   .. ..- attr(*, "lengthvar")= int 3
    ##   .. ..- attr(*, "class")= chr [1:2] "deSolve" "matrix"
    ##   .. ..- attr(*, "type")= chr "rk"
    ##   ..@ initfunc : NULL
