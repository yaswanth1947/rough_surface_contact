rm(list=ls())                ## delete all previous objects
library(RandomFields)
library(spatstat)
library(gstat)
library(sp)
library(spectral)
library(data.table)
library(bit64)
library(reticulate)
library(fractaldim)

#Code to generate random surfaces with Cauchy and Dagum correlation functions. The output file contains the random surface data.

#the output "gen" contains all the RFs

#input
M <- N <- 2048 #grid dimensions
Ns <- 1 #number of realizations 
ext <- 2 #number of extensions for CE
meanv <- 0 #mean value
sigma <- 1 

al <- 1 #alpha in Cauchy correlation function 
bet <- 1.8 #beta in Cauchy correlation function 
#al <- 1 #gamma in Dagum correlation function 
#bet <- 1.8 #beta in Dagum correlation function 

res <- M
step <- 1
W <- list(xrange=c(-((N-1)/2)*step,((N-1)/2)*step),yrange=c(-((N-1)/2)*step,((N-1)/2)*step))

r_gc <- function(u, alpha, beta) (1+abs(u)**alpha)**(-beta/alpha) #Cauchy correlation function'
#r_gc <- function(u, gamma, beta) 1-(1+abs(u)**(-beta))**(-gamma/beta) #Dagum correlation function

grid.prep <- function(W, M, N, ext) {
  cell.width <- diff(W$xrange)/M
  cell.height <- diff(W$yrange)/N
  
  mgrid <- seq(W$xrange[1], W$xrange[2], by = cell.width)
  ngrid <- seq(W$yrange[1], W$yrange[2], by = cell.height)
  mcens <- (mgrid + 0.5 * cell.width)[-(M + 1)]
  ncens <- (ngrid + 0.5 * cell.height)[-(N + 1)]
  
  if (ext <= 1) 
    mgrid.ext <- ngrid.ext <- mcens.ext <- ncens.ext <- M.ext <- N.ext <- NULL else {
      M.ext <- ext * M
      N.ext <- ext * N
      mgrid.ext <- seq(W$xrange[1], W$xrange[2] + (ext - 1) * diff(W$xrange), by = cell.width)
      ngrid.ext <- seq(W$yrange[1], W$yrange[2] + (ext - 1) * diff(W$yrange), by = cell.height)
      mcens.ext <- (mgrid.ext + 0.5 * cell.width)[-(M.ext + 1)]
      ncens.ext <- (ngrid.ext + 0.5 * cell.height)[-(N.ext + 1)]
    }
  
  return(list(M = M, N = N, mgrid = mgrid, ngrid = ngrid, mcens = mcens, ncens = ncens, 
              cell.width = cell.width, cell.height = cell.height, M.ext = M.ext, N.ext = N.ext, 
              mgrid.ext = mgrid.ext, ngrid.ext = ngrid.ext, mcens.ext = mcens.ext, ncens.ext = ncens.ext))
}

timings <- function(ta, tb, tc) {
  tm <- data.frame(tb - ta, tc - tb, tc - ta, row.names = "")
  names(tm) <- c("Decomp", "Generate", "Total elapsed")
  tm
}

covariance.prep <- function(gp,variance,corr.func,...){
  Rx <- gp$M.ext*gp$cell.width
  Ry <- gp$N.ext*gp$cell.height
  m.abs.diff.row1 <- abs(gp$mcens.ext[1]-gp$mcens.ext)
  m.diff.row1 <- pmin(m.abs.diff.row1,Rx-m.abs.diff.row1)
  n.abs.diff.row1 <- abs(gp$ncens.ext[1]-gp$ncens.ext)
  n.diff.row1 <- pmin(n.abs.diff.row1,Ry-n.abs.diff.row1)
  cent.ext.row1 <- expand.grid(m.diff.row1,n.diff.row1)
  D.ext.row1 <- matrix(sqrt(cent.ext.row1[,1]^2+cent.ext.row1[,2]^2),gp$M.ext,gp$N.ext)
  #C.tilde <- variance*r(D.ext.row1,...)
  C.tilde <- variance*r_gc(D.ext.row1,...)
  return(C.tilde)
}

generate.normals <- function(W,M,N,mygrid,Ns,invmat,meanv){
  #realisations <- matrix(NA,prod(dim(invmat)),Ns)
  realisations <- matrix(NA,prod(M,N),Ns)
  cent <- expand.grid(mygrid$mcens, mygrid$ncens)
  for(i in 1:Ns){
    print(i)
    set.seed(i+100)
    realisations[,i] <- -1
    jj <- 0
    field <- invmat*(fft(matrix(rnorm(prod(dim(invmat))),nrow(invmat),ncol(invmat)))/sqrt(prod(dim(invmat))))#randomness here
    realz <- as.vector(Re(fft(field,T)/sqrt(prod(dim(invmat))))[1:M, 1:N])
    realz[!inside.owin(x = cent[, 1], y = cent[, 2], w = W)] <- NA
    realisations[,i] <- matrix(realz, M, N, byrow = TRUE)+meanv
  }
  return(realisations)
}

t1.fft <- t2.1.fft <- c()
for(i in 1:length(res)){
  print(res[i])
  mygrid <- grid.prep(W,res[i],res[i],ext)
  var <- sigma**2
  covmat2 <- covariance.prep(mygrid,var,r_gc,al,bet)
  
  t1.fft[i] <- system.time(
    imat <- sqrt(Re(fft(covmat2,T)))
  )[3]
  
  t2.1.fft[i] <- system.time(
    gen <- generate.normals(W,M,N,mygrid,Ns,imat,meanv)
  )[3]
}

x.seq <- seq(-((N-1)/2)*step, ((N-1)/2)*step, step) 
y.seq <- seq(-((N-1)/2)*step, ((N-1)/2)*step, step) 

aa=expand.grid(x.seq,y.seq)

coords=cbind(aa[,1],aa[,2])
#dataset = matrix(0, nrow = N**2, ncol = Ns+2)
#dataset[,1:2] = coords

#dataset[,3:(Ns+2)] = gen
dataset = matrix(0, nrow = N**2, ncol = Ns)
dataset[,1:(Ns)] = gen

newstring <- paste("cauchyaa",as.character(N),"_",as.character(Ns),"_al",as.character(al),"_b",as.character(bet),'_',as.character(meanv),"cefn",".csv",sep = "")
#newstring <- paste("dagumaa",as.character(N),"_",as.character(Ns),"_al",as.character(al),"_b",as.character(bet),'_',as.character(meanv),"cefn",".csv",sep = "")
#write.csv(array_reshape(dataset, c(M, N)),file=newstring)
write.table(array_reshape(dataset, c(M, N)),file=newstring,sep=",",  col.names=FALSE,  row.names=FALSE)
