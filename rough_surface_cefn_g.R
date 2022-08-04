################### my code ###############################################
rm(list=ls())                ## delete all previous objects
library(RandomFields)        #load the packages before (!) running the script
library(pracma)
library(PracTools) 
library(corTools)

al = 1
bet = 1
model <- RMgencauchy(alpha=al, beta=bet)
#model <- RMdagum(beta=bet, gamma=al)

N = 2048
Del = 0.082 

step <- 1/(N-1) 
x.seq <- seq(0, 1, step) 
y.seq <- seq(0, 1, step) 

for (ii in 1:10){
  print(ii)
  RFoptions(seed=as.integer(ii)) 
  fm <- RFsimulate(model, x = x.seq, y = y.seq, z = NULL, grid=TRUE, spConform=FALSE)
  f <- fm
  newstring <- paste("cauchyaa",as.character(N),"_al",as.character(al),"_b",as.character(bet),"_",num2str(ii),".csv",sep = "")
  write.table(f,file=newstring,sep=",",  col.names=FALSE,  row.names=FALSE)
}
