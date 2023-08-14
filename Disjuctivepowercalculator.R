#disjuctive power finder

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

meantotalpowergivenrunmulipletheta=function(givencomb,A,n,thetas,sd)
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

Disjuctivepowermulipletheta=function(r,n,WA,L,U,sd,thetas)
{
  A=dim(r)[1]
  S=dim(r)[2]
  Allcombs=permutations(n = S,r=A,repeats.allowed = TRUE)
  LAllcombs=dim(Allcombs)[1]
  Dpowerbits=rep(NA,LAllcombs)
  
  CM=correlationmatrixFWER(r,WA) #total correlation matrix
  
  for(iDpower in 1:LAllcombs)
  {
    givencomb=Allcombs[iDpower,] #given combintaion for this run
    #print("givencomb")
    #print(givencomb)
    
    Boundsgivenrun=FWERboundsgivenrun(L,U,givencomb,A)
    Ugiven=Boundsgivenrun$boundU
    Lgiven=Boundsgivenrun$boundL
    #print(Boundsgivenrun)
    #print(CM)
    CMgivenrun=correlationmatrixFWERgivenrun(givencomb,A,S,CM)
    #print(CMgivenrun)
    #should add if dim is 1 as well for 1 arm examples
    
    meangivenrun=meantotalpowergivenrunmulipletheta(givencomb,A,n,thetas,sd)
    #print("meangivenrun")
    #print(meangivenrun)
    if(dim(CMgivenrun)[1]==1)
    {
      set.seed(1)
      Dpowerbits[iDpower]=pmvnorm(lower=Lgiven , upper = Ugiven  ,mean = meangivenrun, sigma = 1)
      
    }
    else
    {
      set.seed(1)
      Dpowerbits[iDpower]=pmvnorm(lower=Lgiven , upper = Ugiven  ,mean = meangivenrun, corr = CMgivenrun)
    }
  }
  return(1-sum(Dpowerbits))
}
load("halfwaydesignPWP.RData")

n=halfwaydesignPWP$n
r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
U=halfwaydesignPWP$U
L=halfwaydesignPWP$L
WA=c(0,n[1])

thetas=c(-log(0.69),-9*10^100)
sd=1
set.seed(1)
PWPsetdis_theta_minf=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)

load("halfwaydesignPWP.RData")

n=halfwaydesignPWP$n
r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
U=halfwaydesignPWP$U
L=halfwaydesignPWP$L
WA=c(0,n[1])

thetas=c(-0,-0)
sd=1
set.seed(1)
PWPsetdis_0_0=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)


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
DPPWL1=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPPWL1

#Line 2
thetas=c(-log(0.69),0)
sd=1
set.seed(1)
DPPWL2=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPPWL2

#Line 3
thetas=c(-log(0.69),-9*10^100)
sd=1
set.seed(1)
DPPWL3=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPPWL3


#Line 4
thetas=c(0,-log(0.69))
sd=1
set.seed(1)
DPPWL4=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPPWL4


#Line 5
thetas=c(0,0)
sd=1
set.seed(1)
DPPWL5=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPPWL5


#Line 6
thetas=c(-9*10^100,-log(0.69))
sd=1
set.seed(1)
DPPWL6=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPPWL6

#conjuctive power
load("halfwaydesignTP.RData")
n=halfwaydesignTP$n
r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
U=halfwaydesignTP$U
L=halfwaydesignTP$L
WA=c(0,n[1])

thetas=c(-log(0.69),-log(0.69))
sd=1
set.seed(1)
DPCPL1=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPCPL1

#Line 2
thetas=c(-log(0.69),0)
sd=1
set.seed(1)
DPCPL2=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPCPL2

#Line 3
thetas=c(-log(0.69),-9*10^100)
sd=1
set.seed(1)
DPCPL3=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPCPL3


#Line 4
thetas=c(0,-log(0.69))
sd=1
set.seed(1)
DPCPL4=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPCPL4


#Line 5
thetas=c(0,0)
sd=1
set.seed(1)
DPCPL5=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPCPL5


#Line 6
thetas=c(-9*10^100,-log(0.69))
sd=1
set.seed(1)
DPCPL6=Disjuctivepowermulipletheta(r,n,WA,L,U,sd,thetas)
DPCPL6
