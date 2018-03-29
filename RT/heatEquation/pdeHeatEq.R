## https://github.com/cwachauf/PDE_solving_R
##
## install.packages("ReacTran")
##
## Examples taken from the book "Solving Differential Equations in R" by Karline Soetart et al.

PDE_Heat_Equation <- function()
{
  require(ReacTran)
  require(deSolve)
  N <- 100
  xgrid <- setup.grid.1D(x.up=0,x.down=1,N=N)
  x <- xgrid$x.mid
  D.coeff <- 0.01
  Diffusion <- function(t,Y,parms)
  {
    tran <- tran.1D(C=Y,C.up=0,C.down=1,D=D.coeff,dx=xgrid)
    list(dY=tran$dC,flux.up=tran$flux.up,flux.down=tran$flux.down)
  }
  ## initial condition:
  Yini <- sin(pi*x)
  times <- seq(from=0,to=5,by=0.01)
  print(system.time(out <- ode.1D(y=Yini,times=times,func=Diffusion,parms=NULL,dimens=N)))
  
  ## plot results:
  par(mfrow=c(1,2))
  plot(out[1,2:(N+1)],x,type="l",lwd=2,xlab="Variable, Y",ylab="Distance, x")
  for(i in seq(2,length(times),by=50))
    lines(out[i,2:(N+1)],x)
  image(out,grid=x,mfrow=NULL,ylab="Distance, x",main="Y")
}

PDE_Diffusion_Advection_Equation <- function()
{
  require(ReacTran)
  require(deSolve)
  N <- 500
  xgrid <- setup.grid.1D(x.up=0,x.down=1,N=N)
  x <- xgrid$x.mid
  D.coeff <- 0.001
  v.coeff <- 0.04
  Diffusion <- function(t,Y,parms)
  {
    tran <- tran.1D(C=Y,C.up=0,C.down=0,D=D.coeff,v=v.coeff,dx=xgrid)
    list(dY=tran$dC,flux.up=tran$flux.up,flux.down=tran$flux.down)
  }
  ## initial condition:
  x0 = 0.2
  ##Yini <-  sin(pi*x)
  Yini <- dnorm(x,mean=0.3,sd=0.05)
  times <- seq(from=0,to=5,by=0.01)
  print(system.time(out <- ode.1D(y=Yini,times=times,func=Diffusion,parms=NULL,dimens=N)))
  
  ## plot results:
  par(mfrow=c(1,2))
  plot(out[1,2:(N+1)],x,type="l",lwd=2,xlab="Variable, Y",ylab="Distance, x")
  for(i in seq(2,length(times),by=50))
    lines(out[i,2:(N+1)],x)
  image(out,grid=x,mfrow=NULL,ylab="Distance, x",main="Y")
}

PDE_Diffusion_Advection_Reaction_Equation <- function()
{
  require(ReacTran)
  require(deSolve)
  N <- 2000
  xgrid <- setup.grid.1D(x.up=0,x.down=5,N=N)
  x <- xgrid$x.mid
  
  ## specify diffusion and advection constants
  DA.coeff <- 0.001
  vA.coeff <- 0.85
  
  DB.coeff <- 0.001
  vB.coeff <- 0.85
  
  DAB.coeff <- 0.0008
  vAB.coeff <- 0.5
  
  k_on <- 0.1
  k_off <- 0.55
  
  ## initial conditions,concentrations:
  cA_0 <- 5
  cB_0 <- 5
  cAB_0 <- 2
  
  Yini_A <- cA_0*dnorm(x,mean=0.2,sd=0.01)
  Yini_B <- cB_0*dnorm(x,mean=0.2,sd=0.01)
  Yini_AB <- cAB_0*dnorm(x,mean=0.2,sd=0.01)
  
  ##Yini <-  sin(pi*x)
  Yini <- c(Yini_A,Yini_B,Yini_AB)
  
  ## times <- seq(from=0,to=7,by=0.01)
  ##print(system.time(out <- ode.1D(y=Yini,times=times,func=Diffusion,parms=NULL,dimens=N)))
  DAR_System <- function(t,Y,parms)
  {
    A <- Y[1:N]
    B <- Y[(N+1):(2*N)]
    AB <- Y[(2*N+1):(3*N)]
    
    dA <- tran.1D(C=A,C.up=0,C.down=0,D=DA.coeff,v=vA.coeff,dx=xgrid)$dC - k_on*A*B + k_off*AB 
    dB <- tran.1D(C=B,C.up=0,C.down=0,D=DB.coeff,v=vB.coeff,dx=xgrid)$dC - k_on*A*B + k_off*AB
    dAB <- tran.1D(C=AB,C.up=0,C.down=0,D=DAB.coeff,v=vAB.coeff,dx=xgrid)$dC + k_on*A*B - k_off*AB
    #dA <- tran.1D(C=A,C.up=0,C.down=0,D=DA.coeff,dx=xgrid)$dC - k_on*A*B + k_off*AB 
    #dB <- tran.1D(C=B,C.up=0,C.down=0,D=DB.coeff,dx=xgrid)$dC - k_on*A*B + k_off*AB
    #dAB <- tran.1D(C=AB,C.up=0,C.down=0,D=DAB.coeff,dx=xgrid)$dC + k_on*A*B - k_off*AB
    list(c(dA,dB,dAB))
    ##tran <- tran.1D(C=Y,C.up=0,C.down=0,D=D.coeff,v=v.coeff,dx=xgrid)
    ##list(dY=tran$dC,flux.up=tran$flux.up,flux.down=tran$flux.down)
    
  }
  ## integrate the stuff...
  times <- seq(from=0,to=3,by=0.001)
  print(system.time(out <- ode.1D(y=Yini,times=times,func=DAR_System,parms=NULL,nspec=3,names=c("A","B","AB"),dimens=N))) 
  image(out,grid=xgrid$x.mid,which="AB")
  return(out)
  ## plot results:
  #par(mfrow=c(1,2))
  #plot(out[1,2:(N+1)],x,type="l",lwd=2,xlab="Variable, Y",ylab="Distance, x")
  #for(i in seq(2,length(times),by=50))
  #  lines(out[i,2:(N+1)],x)
  #image(out,grid=x,mfrow=NULL,ylab="Distance, x",main="Y")
}

Plot_Conc_Profile_At_Index <- function(out,index,npnts)
{
  ymax <- max(out[index,1:(3*npnts)])
  indices <- seq(from=1,to=npnts,by=1)
  plot(indices,out[index,1:npnts],type="l",ylim=c(0,ymax),xlim=c(100,npnts))
  points(indices,out[index,(npnts+1):(2*npnts)],type="l",col="red")
  points(indices,out[index,(2*npnts+1):(3*npnts)],type="l",col="blue")
  
}

Plot_Image_All <- function(out,ntimes,npnts)
{
  out_curr <- out[1:ntimes,1:npnts] + out[1:ntimes,(npnts+1):(2*npnts)] + out[1:ntimes,(2*npnts+1):(3*npnts)]
  image(out_curr)
}