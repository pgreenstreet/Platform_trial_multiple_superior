#pairwisepower
library("mvtnorm")
library("gtools")
kequalkstarPWpower=function(r,gAandS1,gAandS2) #the given AandS1 against the given AandS2
{
  rkj=r[1,gAandS1[2]]
  rkjstar=r[1,gAandS2[2]]
  
  #print("rkj")
  #print(rkj)

  rk=0 #not needed anymoreWA[gAandS1[1]]
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

correlationmatrixPW=function(r)
{
  A=dim(r)[1]
  S=dim(r)[2]
  rowandcol=prod(dim(r)) #how large we need
  Cormatrix=matrix(NA,nrow = S,ncol=S)
  
  AandS1=AandS2=cbind(1,rep(1:S)) #Arm and Stage matrix
  
  for(i in 1:S) #i does each row
  {
    for(istar in 1:S) #istar does each column
    {
      gAandS1=AandS1[i,]
      #print("gAandS1")
      #print(gAandS1)
      gAandS2=AandS2[istar,]
      #print(AandS2)
      #print(gAandS2)
      if(gAandS1[1]==gAandS2[1])
      {
        Cormatrix[i,istar]=kequalkstarPWpower(r,gAandS1,gAandS2)
        #this is always true
      }
      
    }
    
  } 
  
  return(Cormatrix)
}


boundsforPWpower=function(U,L,S)
{
  boundL=vector(mode="list",length = S)
  boundU=vector(mode="list",length = S)
  

  for(i in 1:S)
  {
    boundL[[i]]=c(L[0:(i-1)],U[(i)])
    boundU[[i]]=c(U[0:(i-1)],Inf)
  }
  
  
  return(list("boundL"=boundL,"boundU"=boundU))
}

meanforPWpower=function(n,theta,sd)
{
  means=sqrt(n)*theta/(sd*sqrt(2))
  return(means)
}

pairwisepowereveryarm=function(r,n,theta,sd,L,U)
{
  A=dim(r)[1]
  S=dim(r)[2]
  PWpowerbits=rep(NA,S)
  CM=correlationmatrixPW(r)
  meangivens=meanforPWpower(n,theta,sd)
  Boundsgivenruns=boundsforPWpower(U,L,S)
  #print("Boundsgivenruns")
  #print(Boundsgivenruns)
  Ugivens=Boundsgivenruns$boundU
  Lgivens=Boundsgivenruns$boundL
  for(i in 1:S)
  {
    Uneeded=Ugivens[[i]]
    Lneeded=Lgivens[[i]]
    #print("Lneeded")
    #print(Lneeded)
    CMneeded=CM[1:i,1:i]
    meanneeded=meangivens[1:i]
    if(i==1)
    {
      set.seed(1)
      PWpowerbits[i]=pmvnorm(lower=Lneeded , upper = Uneeded  ,mean = meanneeded, sigma = CMneeded)
      
    }
    else
    {
      set.seed(1)
      PWpowerbits[i]=pmvnorm(lower=Lneeded , upper = Uneeded  ,mean = meanneeded, corr = CMneeded)
    }
  }
  PWpower=sum(PWpowerbits)
  return(PWpower)
}

load("halfwaydesignPWP.RData")
n=halfwaydesignPWP$n
r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
U=halfwaydesignPWP$U
L=halfwaydesignPWP$L
WA=c(0,n[1])

thetas=c(-log(0.69),-log(0.69))
sd=1
set.seed(1)
PWPWL1A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWPWL1A1
set.seed(1)
PWPWL1A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWPWL1A2
#Line 2
thetas=c(-log(0.69),0)
sd=1
set.seed(1)
PWPWL2A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWPWL2A1
set.seed(1)
PWPWL2A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWPWL2A2

#Line 3
thetas=c(-log(0.69),-9*10^100)
sd=1
set.seed(1)
PWPWL3A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWPWL3A1
set.seed(1)
PWPWL3A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWPWL3A2

#Line 4
thetas=c(0,-log(0.69))
sd=1
set.seed(1)
PWPWL4A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWPWL4A1
set.seed(1)
PWPWL4A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWPWL4A2

#Line 5
thetas=c(0,0)
sd=1
set.seed(1)
PWPWL5A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWPWL5A1
set.seed(1)
PWPWL5A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWPWL5A2

#Line 6
thetas=c(-9*10^100,-log(0.69))
sd=1
set.seed(1)
PWPWL6A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWPWL6A1
set.seed(1)
PWPWL6A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWPWL6A2


#now for conjunctive power

load("halfwaydesignTP.RData")
n=halfwaydesignTP$n
r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
U=halfwaydesignTP$U
L=halfwaydesignTP$L
WA=c(0,n[1])

thetas=c(-log(0.69),-log(0.69))
sd=1
set.seed(1)
PWCPL1A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWCPL1A1
set.seed(1)
PWCPL1A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWCPL1A2
#Line 2
thetas=c(-log(0.69),0)
sd=1
set.seed(1)
PWCPL2A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWCPL2A1
set.seed(1)
PWCPL2A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWCPL2A2

#Line 3
thetas=c(-log(0.69),-9*10^100)
sd=1
set.seed(1)
PWCPL3A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWCPL3A1
set.seed(1)
PWCPL3A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWCPL3A2

#Line 4
thetas=c(0,-log(0.69))
sd=1
set.seed(1)
PWCPL4A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWCPL4A1
set.seed(1)
PWCPL4A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWCPL4A2

#Line 5
thetas=c(0,0)
sd=1
set.seed(1)
PWCPL5A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWCPL5A1
set.seed(1)
PWCPL5A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWCPL5A2

#Line 6
thetas=c(-9*10^100,-log(0.69))
sd=1
set.seed(1)
PWCPL6A1=pairwisepowereveryarm(r,n,theta=thetas[1],sd,L,U)
PWCPL6A1
set.seed(1)
PWCPL6A2=pairwisepowereveryarm(r,n,theta=thetas[2],sd,L,U)
PWCPL6A2
