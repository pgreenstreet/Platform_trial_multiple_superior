# Creating FWER for new model... #1 to 1 allocation ratio
library("mvtnorm")
library("gtools")


kequalkstar=function(r,WA,gAandS1,gAandS2) #the given AandS1 against the given AandS2
{
  rkj=r[gAandS1[1],gAandS1[2]]
  rkjstar=r[gAandS2[1],gAandS2[2]]
  
  #print("rkj")
  #print(rkj)
  
  rk=WA[gAandS1[1]]
  r0kj=rkj+rk
  r0kjstar=rkjstar+rk
  
  
  
  part1= (             
          sqrt((rkj)^-1 + (r0kj-rk)^-1)*sqrt((rkjstar)^-1 + (r0kjstar-rk)^-1)    
          )^-1
  
  
  part2op1= ( (rkj)^-1 + (r0kj-rk)^-1) #assuming jstar<j
  
  part2op2= ( (rkjstar)^-1 + (r0kjstar-rk)^-1) #assuming j<jstar
  
  #print("part1")
  #print(part1)
  
  ans=part1*min(part2op1,part2op2)
  #print("ans")
  #print(ans)
  
  return(ans)
}
  

knotequalkstar=function(r,WA,gAandS1,gAandS2) #the given AandS1 against the given AandS2
{
  rkj=r[gAandS1[1],gAandS1[2]]
  rkstarjstar=r[gAandS2[1],gAandS2[2]]
  
  #print("rkj")
  #print(rkj)
  
  rk=WA[gAandS1[1]]
  rkstar=WA[gAandS2[1]]
  r0kj=rkj+rk
  r0kstarjstar=rkstarjstar+rkstar
  
  
  
  part1= (             
    sqrt((rkj)^-1 + (r0kj-rk)^-1)*sqrt((rkstarjstar)^-1 + (r0kstarjstar-rkstar)^-1)    
  )^-1
  
  
  part2= min( (r0kj- max(rk,rkstar)),(r0kstarjstar- max(rk,rkstar)) ) #assuming jstar<j
  
  
  part3=(r0kj-rk)*(r0kstarjstar-rkstar)
  
  #print("part1")
  #print(part1)
  
  ans=max(0,part1*part2/part3)
  #print("ans")
  #print(ans)
  
  return(ans)
}


correlationmatrixFWER=function(r,WA)
{
  A=dim(r)[1]
  S=dim(r)[2]
  rowandcol=prod(dim(r)) #how large we need
  Cormatrix=matrix(NA,nrow = rowandcol,ncol=rowandcol)

  AandS1=AandS2=cbind(rep(1:A,each=S),rep(1:S,times=A)) #Arm and Stage matrix
  
  for(i in 1:rowandcol) #i does each row
  {
    for(istar in 1:rowandcol) #istar does each column
    {
      gAandS1=AandS1[i,]
      #print("gAandS1")
      #print(gAandS1)
      gAandS2=AandS2[istar,]
      if(gAandS1[1]==gAandS2[1])
      {
        Cormatrix[i,istar]=kequalkstar(r,WA,gAandS1,gAandS2)
      }
      else
      {
        Cormatrix[i,istar]=knotequalkstar(r,WA,gAandS1,gAandS2)
      }
      
    }
     
  } 
  
  Cormatrix
  
}


FWERboundsgivenrun=function(L,U,givencomb,A)
{
  boundL=c()
  boundU=c()
  
  for(i in 1:A)
  {
    givencombA=givencomb[i]
    boundLgivenA=c(L[0:(givencombA-1)],-Inf)
    boundUgivenA=c(U[0:(givencombA-1)],L[givencombA])
    boundL=c(boundL,boundLgivenA)
    boundU=c(boundU,boundUgivenA)
  }
  
  return(list("boundL"=boundL,"boundU"=boundU))
  
}

correlationmatrixFWERgivenrun=function(givencomb,A,S,CM)
{
  wneeded=matrix(0,nrow=S,ncol=A) #which are needed
  for(i in 1:A)
  {
    givencombA=givencomb[i]
    wneeded[,i]=c(rep(T,givencombA),rep(F,S-givencombA))
  }
  #print(wneeded)
  #print("c(wneeded)")
  #print(c(wneeded))
  CMrun=CM[c(which(wneeded==1)),c(which(wneeded==1))]
  return(CMrun)
}





#library(MAMS)
#testing against MAMS
#MAMStouse3=mams(K=4, J=3, alpha=0.025, power=0.9, r=1:3, r0=1:3, p=NULL, p0=NULL,ushape="triangular", lshape="triangular",delta=sqrt(2)*qnorm(0.65), delta0=sqrt(2)*qnorm(0.55), sd=1,sample.size = F)
#MAMStouse3


#ans=Boundfinder(r,WA=c(0,0,0,0),amax=20,alpha=0.025,tol = 0.00001)
#ans


#Total power
#the same correlation matrix as before. 
#the means the bit that needs to change.

meanconjuctivepowergivenrun=function(givencomb,A,n,thetas,sd)
{
  
  meangiven=c()
  
  for(i in 1:A)
  {
    givencombA=givencomb[i]
    meangivenA=sqrt(n[1:(givencombA)])*thetas[i]/(sd*sqrt(2))
    
    meangiven=c(meangiven,meangivenA)
  }
  
  return(meangiven)
  
}


conjuctivepowerboundsgivenrun=function(L,U,givencomb,A)
{
  boundL=c()
  boundU=c()
  
  for(i in 1:A)
  {
    givencombA=givencomb[i]
    boundLgivenA=c(L[0:(givencombA-1)],U[givencombA])
    boundUgivenA=c(U[0:(givencombA-1)],Inf)
    boundL=c(boundL,boundLgivenA)
    boundU=c(boundU,boundUgivenA)
  }
  
  return(list("boundL"=boundL,"boundU"=boundU))
  
}

conjuctivepower=function(r,n,thetas,sd,WA,L,U)
{
  A=dim(r)[1]
  S=dim(r)[2]
  Allcombs=permutations(n = S,r=A,repeats.allowed = TRUE)
  LAllcombs=dim(Allcombs)[1]
  conjuctivepowerbits=rep(NA,LAllcombs)
  
  CM=correlationmatrixFWER(r,WA) #total correlation matrix
  
  for(iconjuctivepower in 1:LAllcombs)
  {
    givencomb=Allcombs[iconjuctivepower,] #given combintaion for this run
    #print("givencomb")
    #print(givencomb)
    
    Boundsgivenrun=conjuctivepowerboundsgivenrun(L,U,givencomb,A) #need to change this
    Ugiven=Boundsgivenrun$boundU
    Lgiven=Boundsgivenrun$boundL
    #print(Boundsgivenrun)
    #print(CM)
    CMgivenrun=correlationmatrixFWERgivenrun(givencomb,A,S,CM)
    #print("CMgivenrun")
    #print(CMgivenrun)
    meangivenrun=meanconjuctivepowergivenrun(givencomb,A,n,thetas,sd)
    meangivenrun
    #print("meangivenrun")
    #print(meangivenrun)
    #should add if dim is 1 as well for 1 arm examples
    if(dim(CMgivenrun)[1]==1)
    {
      set.seed(1)
      conjuctivepowerbits[iconjuctivepower]=pmvnorm(lower=Lgiven , upper = Ugiven  ,mean = meangivenrun, sigma = 1)
      
    }
    else
    {
      set.seed(1)
      conjuctivepowerbits[iconjuctivepower]=pmvnorm(lower=Lgiven , upper = Ugiven  ,mean = meangivenrun, corr = CMgivenrun)
    }
  }
  return(sum(conjuctivepowerbits))
}




#pairwisepower
load("halfwaydesignPWP.RData")
n=halfwaydesignPWP$n
r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
U=halfwaydesignPWP$U
L=halfwaydesignPWP$L
WA=c(0,n[1])

thetas=c(-log(0.69),-log(0.69))
sd=1
set.seed(1)
TPPWL1=conjuctivepower(r,n,thetas,sd,WA,L,U)
TPPWL1


#conjunctive power
load("halfwaydesignTP.RData")
n=halfwaydesignTP$n
r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
U=halfwaydesignTP$U
L=halfwaydesignTP$L
WA=c(0,n[1])

thetas=c(-log(0.69),-log(0.69))
sd=1
set.seed(1)
TPCPL1=conjuctivepower(r,n,thetas,sd,WA,L,U)
TPCPL1

