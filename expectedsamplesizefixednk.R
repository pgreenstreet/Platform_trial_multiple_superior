#expected sample size
library("mvtnorm")

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

#FWER


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





meanexpectedsamplesizegivenrun=function(givencomb,A,n,thetaforeach,sd)
{
  
  meangiven=c()
  
  for(i in 1:A)
  {
    givencombA=givencomb[i]
    meangivenA=sqrt(n[1:(givencombA)])*thetaforeach[i]/(sd*sqrt(2))
    
    meangiven=c(meangiven,meangivenA)
  }
  
  return(meangiven)
  
}


totalESSboundsgivenrun=function(L,U,givencomb,givenq,A)
{
  boundL=c()
  boundU=c()
  
  for(i in 1:A)
  {
    givencombA=givencomb[i]
    if(givenq[i]==0)
    {
      LastL=-Inf
      LastU=L[givencombA]
    }
    if(givenq[i]==1)
    {
      LastL=U[givencombA]
      LastU=Inf
    }
    
    boundLgivenA=c(L[0:(givencombA-1)],LastL)
    boundUgivenA=c(U[0:(givencombA-1)],LastU)
    boundL=c(boundL,boundLgivenA)
    boundU=c(boundU,boundUgivenA)
  }
  
  return(list("boundL"=boundL,"boundU"=boundU))
  
}

Expectedsamplesizeprobbit=function(r,n,thetaforeach,sd,WA,L,U)
{
  A=dim(r)[1]
  S=dim(r)[2]
  Allcombs=permutations(n = S,r=A,repeats.allowed = TRUE)
  LAllcombs=dim(Allcombs)[1]
  Allq=permutations(n = 2,v=0:1,r=A,repeats.allowed = TRUE)
  Lallq=dim(Allq)[1]
  expectedsamplesizebits=rep(NA,LAllcombs*Lallq)
  
  CM=correlationmatrixFWER(r,WA) #total correlation matrix
  counter=0
  #print(LAllcombs*Lallq)
  for(iESS in 1:LAllcombs)
  {
    #print(counter)
    givencomb=Allcombs[iESS,] #given combintaion for this run
    #print("givencomb")
    #print(givencomb)
    
  
    #print(Boundsgivenrun)
    #print(CM)
    CMgivenrun=correlationmatrixFWERgivenrun(givencomb,A,S,CM)
    
    meangivenrun=meanexpectedsamplesizegivenrun(givencomb,A,n,thetaforeach,sd)
    
    #print(meangivenrun)
    
    for(iq in 1:Lallq)
    {
      counter=counter+1
      givenq=Allq[iq,]
      givenbounds=totalESSboundsgivenrun(L,U,givencomb,givenq,A)
      #print(givenbounds)
      Lgiven=givenbounds$boundL
      Ugiven=givenbounds$boundU
      if(dim(CMgivenrun)[1]==1)
      {
        set.seed(1)
        expectedsamplesizebits[counter]=pmvnorm(lower=Lgiven, upper = Ugiven  ,mean = meangivenrun, sigma = 1)
          
      }
      else
      {
        set.seed(1)
        expectedsamplesizebits[counter]=pmvnorm(lower=Lgiven, upper = Ugiven  ,mean = meangivenrun, corr = CMgivenrun)
      }
    }
    
  }
  return(expectedsamplesizebits)
}
#n=c(30,60,90)
#WA=c(0,10,40,60)
#r=matrix(c(1,2,3),nrow = 4,ncol=3,byrow=T)
#U=c(2.7,2.5,2)
#L=c(0,0.3,2)
#thetaforeach=c(80,20,-20,-10)
#test=Expectedsamplesize(r,n,thetaforeach,sd=1,WA,L,U)
#sum(test) #should equal 1 

Expectedsamplesizenumberbit=function(n,r,WA)
{
  A=dim(r)[1]
  S=dim(r)[2]
  Allcombs=permutations(n = S,r=A,repeats.allowed = TRUE)
  LAllcombs=dim(Allcombs)[1]
  Allq=permutations(n = 2,v=0:1,r=A,repeats.allowed = TRUE)
  Lallq=dim(Allq)[1]
  
  Activenbit=matrix(rep(NA,LAllcombs*Lallq*A),ncol=A)
  Controlnbit=rep(NA,LAllcombs*Lallq)
  Totalnbit=rep(NA,LAllcombs*Lallq)
  counter=0
  for(iESS in 1:LAllcombs)
  {
    #print(counter)
    givencomb=Allcombs[iESS,]
    for(iq in 1:Lallq)
    {
      counter=counter+1
      givenActiven=n[givencomb]
      givenControln=max(n[givencomb]+WA)
      givenTotaln=sum(givenActiven)+givenControln
      Activenbit[counter,]=givenActiven
      Controlnbit[counter]=givenControln
      Totalnbit[counter]=givenTotaln
    }
    
  }
  return(list("Activen"=Activenbit,"Controln"=Controlnbit,"Totaln"=Totalnbit))
}
#testnbit=Expectedsamplesizenumberbit(n,r,WA)
#sum(testnbit$Totaln*test)


#expected sample size for the pairwise set of data points
ExpectedPWL1=ExpectedPWL2=ExpectedPWL3=ExpectedPWL4=ExpectedPWL5=ExpectedPWL6=rep(NA,201)
load("ResultsPWPdifadding200.RData")

for(runnumber in 1:201)
{
  sd=1
  n=ResultsPWPdifadding200[[runnumber]]$nnow
  r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
  U=ResultsPWPdifadding200[[runnumber]]$U
  L=ResultsPWPdifadding200[[runnumber]]$L
  WA=ResultsPWPdifadding200[[runnumber]]$WA
  
  set.seed(1)
  Npart=Expectedsamplesizenumberbit(n,r,WA)
  
  #line1
  thetaforeach=c(-log(0.69),-log(0.69))
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedPWL1[runnumber]=sum(Ppart*Npart$Totaln)
  
  
  #line2
  thetaforeach=c(-log(0.69),0)
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedPWL2[runnumber]=sum(Ppart*Npart$Totaln)
  
  #line3
  thetaforeach=c(-log(0.69),-9*10^100)
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedPWL3[runnumber]=sum(Ppart*Npart$Totaln)
  
  #line4
  thetaforeach=c(0,-log(0.69))
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedPWL4[runnumber]=sum(Ppart*Npart$Totaln)
  
  #line5
  thetaforeach=c(0,0)
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedPWL5[runnumber]=sum(Ppart*Npart$Totaln)
  
  #line6
  thetaforeach=c(-9*10^100,-log(0.69))
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedPWL6[runnumber]=sum(Ppart*Npart$Totaln)
  
  print(runnumber)
}
ESSeachrunPWPall200=list("line1"=ExpectedPWL1,"line2"=ExpectedPWL2,"line3"=ExpectedPWL3,
                      "line4"=ExpectedPWL4,"line5"=ExpectedPWL5,"line6"=ExpectedPWL6)
save(ESSeachrunPWPall200,file = "ESSeachrunPWPall200.RData")


#expected sample size for the pairwise set of data points
ExpectedTPL1=ExpectedTPL2=ExpectedTPL3=ExpectedTPL4=ExpectedTPL5=ExpectedTPL6=rep(NA,201)
load("ResultsTPdifadding200.RData")

for(runnumber in 1:201)
{
  sd=1
  n=ResultsTPdifadding200[[runnumber]]$nnow
  r=matrix(c(1,2,1,2),nrow = 2,byrow = T)*n[1]
  U=ResultsTPdifadding200[[runnumber]]$U
  L=ResultsTPdifadding200[[runnumber]]$L
  WA=ResultsTPdifadding200[[runnumber]]$WA
  
  set.seed(1)
  Npart=Expectedsamplesizenumberbit(n,r,WA)
  
  #line1
  thetaforeach=c(-log(0.69),-log(0.69))
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedTPL1[runnumber]=sum(Ppart*Npart$Totaln)
  
  
  #line2
  thetaforeach=c(-log(0.69),0)
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedTPL2[runnumber]=sum(Ppart*Npart$Totaln)
  
  #line3
  thetaforeach=c(-log(0.69),-9*10^100)
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedTPL3[runnumber]=sum(Ppart*Npart$Totaln)
  
  #line4
  thetaforeach=c(0,-log(0.69))
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedTPL4[runnumber]=sum(Ppart*Npart$Totaln)
  
  #line5
  thetaforeach=c(0,0)
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedTPL5[runnumber]=sum(Ppart*Npart$Totaln)
  
  #line6
  thetaforeach=c(-9*10^100,-log(0.69))
  set.seed(1)
  Ppart=Expectedsamplesizeprobbit(r,n,thetaforeach,sd,WA,L,U)
  ExpectedTPL6[runnumber]=sum(Ppart*Npart$Totaln)
  
  print(runnumber)
}
ESSeachrunTPall200=list("line1"=ExpectedTPL1,"line2"=ExpectedTPL2,"line3"=ExpectedTPL3,
                     "line4"=ExpectedTPL4,"line5"=ExpectedTPL5,"line6"=ExpectedTPL6)
save(ESSeachrunTPall200,file = "ESSeachrunTPall200.RData")
