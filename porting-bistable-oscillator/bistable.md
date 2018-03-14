simecol bistable
================

``` r
library("simecol")
```

    ## Loading required package: deSolve

``` r
bis.model <- new("odeModel",                   # object class = 'odeModel'
                
         main = function(times, y, parms) {   # object method
                                              # simecol magic: 2nd argument 'y' points to 'init'
          with(as.list(c(parms, y)), {  
            
            change <- kineticEnergy
            change0 <- factor * ( potentialEnergy - potentialEnergy^3 )  - friction*kineticEnergy
            kineticEnergy <- change0
            potentialEnergy <- change
            list(c(kineticEnergy, potentialEnergy))    # object attributes (data) 
           })
         },
         
        times = seq(0, 50, .015625),
        parms = c(startingPoint=1, friction=0.01, factor=1),
        init = c(kineticEnergy=1, potentialEnergy=1),
        solver = "rk4"                      # 'lsoda' available for 'odemodel' class
        
      )
```

``` r
bis.sim <- sim(bis.model)
plot(bis.sim)
```

![](bistable_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
plotupca <- function(obj, ...) {
  
  library(gridExtra)
  library(ggplot2)
  
  o.df <- out(obj)  # output a dataframe
  
  gg1 <- ggplot()  +
    geom_path(aes(x=time, y=potentialEnergy), o.df,  
        colour = "purple", size=1, alpha=0.6)  +
    geom_path(aes(x=time, y=kineticEnergy), o.df, colour = "green", size=1, alpha=0.8)
  
  gg2 <- ggplot()  +
    geom_path(aes(x=potentialEnergy, y=kineticEnergy), o.df,  colour = "purple", alpha=0.6) 
  
  grid.arrange(gg1, gg2)

}
```

``` r
 plotupca(bis.sim)
```

![](bistable_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
head(out(bis.sim))
```

    ##       time kineticEnergy potentialEnergy
    ## 1 0.000000     1.0000000        1.000000
    ## 2 0.015625     0.9995958        1.015622
    ## 3 0.031250     0.9986806        1.031235
    ## 4 0.046875     0.9972315        1.046828
    ## 5 0.062500     0.9952258        1.062395
    ## 6 0.078125     0.9926407        1.077926

``` r
main(bis.sim)
```

    ## function(times, y, parms) {   # object method
    ##                                               # simecol magic: 2nd argument 'y' points to 'init'
    ##           with(as.list(c(parms, y)), {  
    ##             
    ##             change <- kineticEnergy
    ##             change0 <- factor * ( potentialEnergy - potentialEnergy^3 )  - friction*kineticEnergy
    ##             kineticEnergy <- change0
    ##             potentialEnergy <- change
    ##             list(c(kineticEnergy, potentialEnergy))    # object attributes (data) 
    ##            })
    ##          }

``` r
init(bis.sim)
```

    ##   kineticEnergy potentialEnergy 
    ##               1               1

``` r
parms(bis.sim)
```

    ## startingPoint      friction        factor 
    ##          1.00          0.01          1.00

``` r
equations(bis.sim)
```

    ## NULL

``` r
solver(bis.sim)
```

    ## [1] "rk4"

``` r
class(bis.sim)
```

    ## [1] "odeModel"
    ## attr(,"package")
    ## [1] "simecol"

``` r
str(bis.sim)
```

    ## Formal class 'odeModel' [package "simecol"] with 10 slots
    ##   ..@ parms    : Named num [1:3] 1 0.01 1
    ##   .. ..- attr(*, "names")= chr [1:3] "startingPoint" "friction" "factor"
    ##   ..@ init     : Named num [1:2] 1 1
    ##   .. ..- attr(*, "names")= chr [1:2] "kineticEnergy" "potentialEnergy"
    ##   ..@ observer : NULL
    ##   ..@ main     :function (times, y, parms)  
    ##   .. ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 6 17 16 10 17 10 6 16
    ##   .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x7f8ef9d9e778> 
    ##   ..@ equations: NULL
    ##   ..@ times    : num [1:3201] 0 0.0156 0.0312 0.0469 0.0625 ...
    ##   ..@ inputs   : NULL
    ##   ..@ solver   : chr "rk4"
    ##   ..@ out      : deSolve [1:3201, 1:3] 0 0.0156 0.0312 0.0469 0.0625 ...
    ##   .. ..- attr(*, "istate")= int [1:21] 0 3200 12801 NA NA NA NA NA NA NA ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   .. .. ..$ : NULL
    ##   .. .. ..$ : chr [1:3] "time" "kineticEnergy" "potentialEnergy"
    ##   .. ..- attr(*, "lengthvar")= int 2
    ##   .. ..- attr(*, "class")= chr [1:2] "deSolve" "matrix"
    ##   .. ..- attr(*, "type")= chr "rk"
    ##   ..@ initfunc : NULL
