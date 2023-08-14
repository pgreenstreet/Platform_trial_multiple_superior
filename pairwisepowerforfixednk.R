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



# now looking at finding n in the case where WA is stating how many after you add the treatment

#FWER functions needed
tri=function(rf) 
{
  rfNA=rf[!is.na(rf)]
  MrfNA=max(rfNA)
  LrfNA=length(rfNA)
  rfNA1=rfNA[1]
  rfNAa1=rfNA/rfNA1
  MrfNAa1=MrfNA/rfNA1
  UfNA=1*(1+rfNAa1/MrfNAa1)/sqrt(rfNAa1)
  LfNA=-1*(1-3*rfNAa1/MrfNAa1)/sqrt(rfNAa1)
  
  NNA=length(rf)-LrfNA #number of NA needed
  Uf=c(UfNA,rep(NA,NNA))
  Lf=c(LfNA,rep(NA,NNA))
  
  return(list("U"=Uf,"L"=Lf))  
}


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


FWER=function(r,WA,L,U)
{
  A=dim(r)[1]
  S=dim(r)[2]
  Allcombs=permutations(n = S,r=A,repeats.allowed = TRUE)
  LAllcombs=dim(Allcombs)[1]
  FWERbits=rep(NA,LAllcombs)
  
  CM=correlationmatrixFWER(r,WA) #total correlation matrix
  
  for(iFWER in 1:LAllcombs)
  {
    givencomb=Allcombs[iFWER,] #given combintaion for this run
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
    if(dim(CMgivenrun)[1]==1)
    {
      set.seed(1)
      FWERbits[iFWER]=pmvnorm(lower=Lgiven , upper = Ugiven  ,mean = rep(0,dim(CMgivenrun)[1]), sigma = 1)
      
    }
    else
    {
      set.seed(1)
      FWERbits[iFWER]=pmvnorm(lower=Lgiven , upper = Ugiven  ,mean = rep(0,dim(CMgivenrun)[1]), corr = CMgivenrun)
    }
  }
  return(1-sum(FWERbits))
}


Boundfinder=function(r,WA,amax,alpha,tol)
{
  bondshapewanted=tri(r[1,]) #need to make these functions and change this for the differnt ones
  
  cFWER=-Inf#current FWER
  ua=2*amax
  la=0
  ca=(ua+la)/2
  while( cFWER>(alpha) | cFWER<(alpha-tol) )
  {
    if(cFWER>(alpha))
    {
      ua=ua
      la=ca
      ca=(ua+la)/2
      cU=ca*bondshapewanted$U
      cL=ca*bondshapewanted$L
      cFWER=FWER(r = r,WA = WA,L = cL,U = cU)
      #print(cFWER)
    }
    
    if(cFWER<(alpha))
    {
      ua=ca
      la=la
      ca=(ua+la)/2
      cU=ca*bondshapewanted$U
      cL=ca*bondshapewanted$L
      cFWER=FWER(r = r,WA = WA,L = cL,U = cU)
      #print(cFWER)
    }
  }
  U=ca*bondshapewanted$U
  L=ca*bondshapewanted$L
  return(list("L"=L,"U"=U))
}

####

simplenfinderPWpower=function(r,theta,sd,amax,nmax,alpha,beta,tol) #finds it when everystarts at ones
{
  WA=rep(0,dim(r)[1])
  simplebounds=Boundfinder(r,WA,amax,alpha,tol)
  L=simplebounds$L
  U=simplebounds$U
  r=r #/GCF(r)#first step is to reduce r to its simplest form #have not bothered with this
  cPWpower=Inf #current PW
  un=2*nmax
  ln=0
  cn=(un+ln)/2
  
  while( cPWpower>(1-beta+tol) | cPWpower<(1-beta) )
  {
    if(cPWpower>(1-beta))
    {
      un=cn
      ln=ln
      cn=(un+ln)/2
      testn=cn*r[1,] #as it the same for every row in the code
      
      cPWpower=pairwisepowereveryarm(r=r,n=testn,theta=theta,sd=sd,L=L,U=U)
      #print(cPWpower)
    }
    if(cPWpower<(1-beta))
    {
      un=un
      ln=cn
      cn=(un+ln)/2
      testn=cn*r[1,] #as it the same for every row in the code
      cPWpower=pairwisepowereveryarm(r=r,n=testn,theta=theta,sd=sd,L=L,U=U)
      #print(cPWpower)
    }
  }
  #print(testn)
  n=ceiling(cn)*r[1,]
  cPWpower=pairwisepowereveryarm(r=r,n=n,theta=theta,sd=sd,L=L,U=U)
  #print(cPWpower)
  return(list("n"=n,"U"=U,"L"=L))
}


nfinderPWpowerwhensetadding=function(r,theta,sd,WA,amax,nmax,alpha,beta,tol)
{
  nnow=simplenfinderPWpower(r=r,theta,sd,amax,nmax,alpha,beta,tol)$n
  nnow=nnow-r[1,]
  
  A=dim(r)[1]
  S=dim(r)[2]
  
  cPWpower=0
  while(cPWpower<(1-beta))
  {
    nnow=nnow+r[1,]
    rnow=matrix(nnow,nrow = A,ncol=S,byrow=T)
    newbounds=Boundfinder(rnow,WA,amax,alpha,tol)
    L=newbounds$L
    U=newbounds$U
    cPWpower=pairwisepowereveryarm(r=rnow,n=nnow,theta=theta,sd=sd,L=L,U=U)
    #print(cPWpower)
  }
  
  return(list("L"=L,"U"=U,"nnow"=nnow,"WA"=WA)) 
}

WA=c(0,0) #States when treatments are added to the trial
r=matrix(c(1,2),nrow = 2,ncol=2,byrow=T) # Change this for the number of arms and stages
nfinderPWpowerwhensetadding(r=r,theta = -log(0.69),sd = 1,WA = WA,amax = 25,nmax = 1000,alpha = 0.025,beta = 0.2,tol = 0.00001)

ResultsPWPdifadding200=vector(mode="list",length = 201)

for(WAsv in 0:200) #WA second value
{
  cWA=c(0,WAsv)
  set.seed(1)
  ResultsPWPdifadding200[[WAsv+1]]=nfinderPWpowerwhensetadding(r=r,theta = -log(0.69),sd = 1,WA = cWA,amax = 25,nmax = 1000,alpha = 0.025,beta = 0.2,tol = 0.00001)
  print(WAsv)
}  

save(ResultsPWPdifadding200,file="ResultsPWPdifadding200.RData")
#maximum sample size finder
maxeachrunPWP200=rep(NA,201)
for(i in 1:201)
{
  maxeachrunPWP200[i]=(ResultsPWPdifadding200[[i]]$nnow[2])*3+ResultsPWPdifadding200[[i]]$WA[2]
}

plot(y=maxeachrunPWP200,x=0:200)
save(maxeachrunPWP200,file="maxeachrunPWP200.RData")
