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
            
            # kineticEnergy = change0
            # potentialEnergy = change

            change0dt <- change0(friction, factor, potentialEnergy, kineticEnergy, times)
            changedt <- kineticEnergy
            
            kineticEnergy <- changedt
            potentialEnergy <- change0dt
            
            list(c(kineticEnergy, potentialEnergy))    # object attributes (data) 
            
            # change0 = potential energy - POWER(potential energy ,3)  - friction*kinetic energy + Excitation
            # change = kinetic energy
            # Excitation = SIN(Time)
            
           })
         },
         equations = list(
           excite = function(i) {sin( i / (2*pi) )},
           power = function(i) {i^3},
           change = function(i) {i},
           change0 = function(h, i, j, k, t) { (j - power(j) ) - h * k + excite(t)}
         ),
        times = seq(0, 50, .015625),
        parms = c(factor=1, friction=.01),
        init = c(kineticEnergy=1, potentialEnergy=1),
        solver = "lsoda"                      # 'lsoda' available for 'odemodel' class
        
      )
```

``` r
bis.sim <- sim(bis.model)
plot(bis.sim)
```

![](bistableOscExc_files/figure-markdown_github/unnamed-chunk-2-1.png)

``` r
plotupca <- function(obj, ...) {
  
  library(gridExtra)
  library(ggplot2)
  
  o.df <- out(obj)  # output a dataframe
  
  gg1 <- ggplot()  +
    geom_line(aes(x=time, y=potentialEnergy), o.df, colour = "purple", size=1, alpha=0.8)  +
    geom_line(aes(x=time, y=kineticEnergy), o.df, colour = "green", size=1, alpha=0.8)
  
  gg2 <- ggplot()  +
    geom_path(aes(x=potentialEnergy, y=kineticEnergy), o.df,  colour = "purple", alpha=0.6) 
  
  
  grid.arrange(gg1, gg2, ncol=2, nrow=1)

}
```

``` r
 plotupca(bis.sim)
```

![](bistableOscExc_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
head(out(bis.sim))
```

    ##       time kineticEnergy potentialEnergy
    ## 1 0.000000      1.000000       1.0000000
    ## 2 0.015625      1.015748       0.9998641
    ## 3 0.031250      1.031744       0.9997681
    ## 4 0.046875      1.047992       0.9997109
    ## 5 0.062500      1.064495       0.9996914
    ## 6 0.078125      1.081259       0.9997081

``` r
main(bis.sim)
```

    ## function(times, y, parms) {   # object method
    ##                                               # simecol magic: 2nd argument 'y' points to 'init'
    ##           with(as.list(c(parms, y)), {
    ##             
    ##             # kineticEnergy = change0
    ##             # potentialEnergy = change
    ## 
    ##             change0dt <- change0(friction, factor, potentialEnergy, kineticEnergy, times)
    ##             changedt <- kineticEnergy
    ##             
    ##             kineticEnergy <- changedt
    ##             potentialEnergy <- change0dt
    ##             
    ##             list(c(kineticEnergy, potentialEnergy))    # object attributes (data) 
    ##             
    ##             # change0 = potential energy - POWER(potential energy ,3)  - friction*kinetic energy + Excitation
    ##             # change = kinetic energy
    ##             # Excitation = SIN(Time)
    ##             
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

    ##   factor friction 
    ##     1.00     0.01

``` r
equations(bis.sim)
```

    ## $excite
    ## function (i) 
    ## {
    ##     sin(i/(2 * pi))
    ## }
    ## 
    ## $power
    ## function (i) 
    ## {
    ##     i^3
    ## }
    ## 
    ## $change
    ## function (i) 
    ## {
    ##     i
    ## }
    ## 
    ## $change0
    ## function (h, i, j, k, t) 
    ## {
    ##     (j - power(j)) - h * k + excite(t)
    ## }

``` r
solver(bis.sim)
```

    ## [1] "lsoda"

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
    ##   ..@ parms    : Named num [1:2] 1 0.01
    ##   .. ..- attr(*, "names")= chr [1:2] "factor" "friction"
    ##   ..@ init     : Named num [1:2] 1 1
    ##   .. ..- attr(*, "names")= chr [1:2] "kineticEnergy" "potentialEnergy"
    ##   ..@ observer : NULL
    ##   ..@ main     :function (times, y, parms)  
    ##   .. ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 6 17 26 10 17 10 6 26
    ##   .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x7fd3f45f3348> 
    ##   ..@ equations:List of 4
    ##   .. ..$ excite :function (i)  
    ##   .. .. ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 28 21 28 51 21 51 28 28
    ##   .. .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x7fd3f45f3348> 
    ##   .. ..$ power  :function (i)  
    ##   .. .. ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 29 20 29 36 20 36 29 29
    ##   .. .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x7fd3f45f3348> 
    ##   .. ..$ change :function (i)  
    ##   .. .. ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 30 21 30 35 21 35 30 30
    ##   .. .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x7fd3f45f3348> 
    ##   .. ..$ change0:function (h, i, j, k, t)  
    ##   .. .. ..- attr(*, "srcref")=Class 'srcref'  atomic [1:8] 31 22 31 83 22 83 31 31
    ##   .. .. .. .. ..- attr(*, "srcfile")=Classes 'srcfilecopy', 'srcfile' <environment: 0x7fd3f45f3348> 
    ##   ..@ times    : num [1:3201] 0 0.0156 0.0312 0.0469 0.0625 ...
    ##   ..@ inputs   : NULL
    ##   ..@ solver   : chr "lsoda"
    ##   ..@ out      : deSolve [1:3201, 1:3] 0 0.0156 0.0312 0.0469 0.0625 ...
    ##   .. ..- attr(*, "dimnames")=List of 2
    ##   .. .. ..$ : NULL
    ##   .. .. ..$ : chr [1:3] "time" "kineticEnergy" "potentialEnergy"
    ##   .. ..- attr(*, "istate")= int [1:21] 2 3470 5038 NA 3 3 0 52 22 NA ...
    ##   .. ..- attr(*, "rstate")= num [1:5] 0.0156 0.0156 50.0035 0 11.6216
    ##   .. ..- attr(*, "lengthvar")= int 2
    ##   .. ..- attr(*, "class")= chr [1:2] "deSolve" "matrix"
    ##   .. ..- attr(*, "type")= chr "lsoda"
    ##   ..@ initfunc : NULL
