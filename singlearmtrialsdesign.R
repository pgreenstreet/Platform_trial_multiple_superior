#simple 1 stage trial.

FWER=function(L,U)
{
  set.seed(1)
  FWERpart1=pmvnorm(lower=U[1], upper = Inf  ,mean = 0, sigma = 1)
  set.seed(1)
  FWERpart2=pmvnorm(lower=c(L[1],U[2]), upper = c(U[1],Inf)  ,mean = c(0,0),
                    corr = matrix(c(1,sqrt(.5),sqrt(.5),1),nrow = 2,byrow=T) )
  FWERfound=FWERpart1+FWERpart2
  return(FWERfound[[1]])
}
#FWER(c(0,2),c(2.4,2.4))


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

Boundfinder=function(amax,alpha,tol)
{
  bondshapewanted=tri(c(1,2)) #need to make these functions and change this for the differnt ones
  
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
      cFWER=FWER(L = cL,U = cU)
      #print(cFWER)
    }
    
    if(cFWER<(alpha))
    {
      ua=ca
      la=la
      ca=(ua+la)/2
      cU=ca*bondshapewanted$U
      cL=ca*bondshapewanted$L
      cFWER=FWER(L = cL,U = cU)
      #print(cFWER)
    }
  }
  U=ca*bondshapewanted$U
  L=ca*bondshapewanted$L
  return(list("L"=L,"U"=U))
}
#testbounds=Boundfinder(10,0.025,0.00001)

powerfinder=function(L,U,delta,n,sd)
{
  meansneeded=delta*sqrt(n)/(sd*sqrt(2))
  set.seed(1)
  Powerpart1=pmvnorm(lower=U[1], upper = Inf  ,mean = meansneeded[1], sigma = 1)
  set.seed(1)
  Powerpart2=pmvnorm(lower=c(L[1],U[2]), upper = c(U[1],Inf)  ,mean = meansneeded,
                    corr = matrix(c(1,sqrt(.5),sqrt(.5),1),nrow = 2,byrow=T) )
  Powerfound=Powerpart1+Powerpart2
  return(Powerfound[[1]])
}

samplesizefounder=function(L,U,delta,beta,sd)
{
  nnow=c(0,0)
  cpower=0
  while(cpower<=(1-beta))
  {
    nnow=nnow+c(1,2)
    cpower=powerfinder(L=L,U=U,delta=delta,n=nnow,sd=sd)
    #print(cpower)
  }
  return(nnow)
}
#samplesizefounder(L=testbounds$L,U=testbounds$U,delta = -log(0.69),beta = 0.2,sd=1)

Onearmdesign=function(alpha,beta,amax,tol,delta,sd)
{
  Thebounds=Boundfinder(amax,alpha,tol)
  Then=samplesizefounder(L=Thebounds$L,U=Thebounds$U,delta,beta,sd)
  return(list("L"=Thebounds$L,"U"=Thebounds$U,"n"=Then))
}


ESSonearm=function(L,U,n,delta,sd)
{
  meansneeded=delta*sqrt(n)/(sd*sqrt(2))
  #first stage stop
  set.seed(1)
  part1stop1=pmvnorm(lower=U[1], upper = Inf  ,mean = meansneeded[1], sigma = 1)
  set.seed(1)
  part2stop1=pmvnorm(lower=-Inf, upper = L[1]  ,mean = meansneeded[1], sigma = 1)
  #print("n")
  #print(n)
  Stage1stop=(part1stop1+part2stop1)*2*n[1]
  
  #print("(part1stop1+part2stop1)")
  #print((part1stop1+part2stop1))
  #second stage stop
  set.seed(1)
  part1stop2=pmvnorm(lower=c(L[1],U[2]), upper = c(U[1],Inf)  ,mean = meansneeded,
                     corr = matrix(c(1,sqrt(.5),sqrt(.5),1),nrow = 2,byrow=T) )
  set.seed(1)
  part2stop2=pmvnorm(lower=c(L[1],-Inf), upper = c(U[1],U[2])  ,mean = meansneeded,
                     corr = matrix(c(1,sqrt(.5),sqrt(.5),1),nrow = 2,byrow=T) )
  #print("(part1stop2+part2stop2)")
  #print((part1stop2+part2stop2))
  Stage2stop=(part1stop2+part2stop2)*2*n[2]
  #print(((part1stop2+part2stop2)+(part1stop1+part2stop1))[[1]])
  TotalESS=Stage1stop[[1]]+Stage2stop[[1]]
  return(TotalESS)
}
Designfor1armPWP=Onearmdesign(0.025,0.2,20,0.00001,-log(0.69),1)
Designfor1armPWP
delta0ESSPWP=ESSonearm(Designfor1armPWP$L,Designfor1armPWP$U,n=Designfor1armPWP$n,delta=0,sd=1)

deltalog.69ESSPWP=ESSonearm(Designfor1armPWP$L,Designfor1armPWP$U,n=Designfor1armPWP$n,delta=-log(0.69),sd=1)

deltalogmInfESSPWP=ESSonearm(Designfor1armPWP$L,Designfor1armPWP$U,n=Designfor1armPWP$n,delta=--9*10^100,sd=1)

MaxsepatePWP=Designfor1armPWP$n[2]*4
Line1sepatePWP=deltalog.69ESSPWP+deltalog.69ESSPWP
Line2sepatePWP=deltalog.69ESSPWP+delta0ESSPWP
Line3sepatePWP=deltalog.69ESSPWP+deltalogmInfESSPWP
Line4sepatePWP=delta0ESSPWP+deltalog.69ESSPWP
Line5sepatePWP=delta0ESSPWP+delta0ESSPWP
Line6sepatePWP=deltalogmInfESSPWP+deltalog.69ESSPWP

SamplesizesepatePWP=list("MaxsepatePWP"=MaxsepatePWP,"Line1sepatePWP"=Line1sepatePWP,"Line2sepatePWP"=Line2sepatePWP,
                        "Line3sepatePWP"=Line3sepatePWP,"Line4sepatePWP"=Line4sepatePWP,"Line5sepatePWP"=Line5sepatePWP
                        ,"Line6sepatePWP"=Line6sepatePWP)
save(SamplesizesepatePWP,file="SamplesizesepatePWP.RData")


Designfor1armTP=Onearmdesign(0.025,1-sqrt(.8),20,0.00001,-log(0.69),1)
Designfor1armTP

delta0ESSTP=ESSonearm(Designfor1armTP$L,Designfor1armTP$U,n=Designfor1armTP$n,delta=0,sd=1)

deltalog.69ESSTP=ESSonearm(Designfor1armTP$L,Designfor1armTP$U,n=Designfor1armTP$n,delta=-log(0.69),sd=1)

deltalogmInfESSTP=ESSonearm(Designfor1armTP$L,Designfor1armTP$U,n=Designfor1armTP$n,delta=--9*10^100,sd=1)
MaxsepateTP=Designfor1armTP$n[2]*4
Line1sepateTP=deltalog.69ESSTP+deltalog.69ESSTP
Line2sepateTP=deltalog.69ESSTP+delta0ESSTP
Line3sepateTP=deltalog.69ESSTP+deltalogmInfESSTP
Line4sepateTP=delta0ESSTP+deltalog.69ESSTP
Line5sepateTP=delta0ESSTP+delta0ESSTP
Line6sepateTP=deltalogmInfESSTP+deltalog.69ESSTP

SamplesizesepateTP=list("MaxsepateTP"=MaxsepateTP,"Line1sepateTP"=Line1sepateTP,"Line2sepateTP"=Line2sepateTP,
                        "Line3sepateTP"=Line3sepateTP,"Line4sepateTP"=Line4sepateTP,"Line5sepateTP"=Line5sepateTP
                        ,"Line6sepateTP"=Line6sepateTP)
save(SamplesizesepateTP,file="SamplesizesepateTP.RData")
