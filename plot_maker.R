#for graph TP i need 196 for PWP I need 154
#creating the figure
#pairwisepower
load(file = "ESSeachrunPWPall200.RData")
load(file="maxeachrunPWP200.RData")
minPW=min(ESSeachrunPWPall200$line5[1:155]) #minimum value
maxPW=max(maxeachrunPWP200[1:155]) #minimum value
length(maxeachrunPWP200[1:155])
#plot
plot(x=0:154,y=maxeachrunPWP200[1:155],xlim=c(0,154),ylim = c(minPW-50,maxPW+10) ,main="Setting 2", type="l",col="Black",lty=1,lwd = 2,xlab = "When the second active treatment is added",ylab = "Required sample size")

#name1=expression("E(N|" (theta[1] == 0)) expression((theta[1] == 0))

name1=expression("max(N)")
name2=expression("E(N|"~(theta~"',"~ theta~"'")~")" ) 
name3=expression("E(N|"~(theta~"',"~ 0)~")" ) 
name4=expression("E(N|"~(theta~"',"~ -infinity)~")" ) 
name5=expression("E(N|"~(0~","~ theta~"'")~")" ) 
name6=expression("E(N|"~(0~","~ 0)~")" ) 
name7=expression("E(N|"~(-infinity ~","~ theta~"'")~")" ) 


legend(100, 370, legend=c(name1, name2, name3, name4,name5,name6,name7 ),
       col=c("Black","Dark Blue","Dark Green","Purple","Orange","Red","Blue"),lwd = 2, cex=0.8,
       box.lty=0)

points(x=0:154,y=ESSeachrunPWPall200$line1[1:155],type="l",col="Dark Blue",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line2[1:155],type="l",col="Dark Green",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line3[1:155],type="l",col="Purple",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line4[1:155],type="l",col="Orange",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line5[1:155],type="l",col="Red",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line6[1:155],type="l",col="Blue",lwd = 2,lty=1)


name2=expression(alpha[sub1^"'"*sub2])

load(file="SamplesizesepatePWP.RData")

SamplesizesepatePWP$MaxsepatePWP
#abline(h=SamplesizesepatePWP$MaxsepatePWP,lty=3,col="black")
#abline(h=SamplesizesepatePWP$Line1sepatePWP,lty=3,col="Dark Blue")
#need to find when this goes above


nopointmax=min(which(SamplesizesepatePWP$MaxsepatePWP<=maxeachrunPWP200))-1
nopointline1=min(which(SamplesizesepatePWP$Line1sepatePWP<=ESSeachrunPWPall200$line1))-1
nopointline2=min(which(SamplesizesepatePWP$Line2sepatePWP<=ESSeachrunPWPall200$line2))-1
nopointline3=min(which(SamplesizesepatePWP$Line3sepatePWP<=ESSeachrunPWPall200$line3))-1
nopointline4=min(which(SamplesizesepatePWP$Line4sepatePWP<=ESSeachrunPWPall200$line4))-1
nopointline5=min(which(SamplesizesepatePWP$Line5sepatePWP<=ESSeachrunPWPall200$line5))-1
nopointline6=min(which(SamplesizesepatePWP$Line6sepatePWP<=ESSeachrunPWPall200$line6))-1


abline(v=nopointmax,lty=3,lwd = 3,col="Black")
abline(v=nopointline1,lty=3,lwd = 3,col="Dark Blue")
abline(v=nopointline2,lty=3,lwd = 3,col="Dark Green")
abline(v=nopointline3,lty=3,lwd = 3,col="Purple")
abline(v=nopointline4,lty=3,lwd = 3,col="Orange")
abline(v=nopointline5,lty=3,lwd = 3,col="Red")
abline(v=nopointline6,lty=3,lwd = 3,col="Blue")
#5.54 by 6.72
#PWPplotS2
#PWPplotSSR2

#values for paper 
min(maxeachrunPWP200[1:155])
max(maxeachrunPWP200[1:155])
nopointmax
SamplesizesepatePWP$MaxsepatePWP

nopointline6
SamplesizesepatePWP$Line6sepatePWP

#totalpower
load(file = "ESSeachrunTPall200.RData")
load(file="maxeachrunTP200.RData")
minTP=min(ESSeachrunTPall200$line5[1:197]) #minimum value
maxTP=max(maxeachrunTP200[1:197]) #minimum value
length(maxeachrunTP200[1:197])
#plot
plot(x=0:196,y=maxeachrunTP200[1:197],xlim=c(0,196),ylim = c(minTP-50,maxTP+10) ,type="l",main="Setting 2",col="Black",lty=1,lwd = 2,xlab = "When the second active treatment is added",ylab = "Required sample size")

#name1=expression("E(N|" (theta[1] == 0)) expression((theta[1] == 0))

name1=expression("max(N)")
name2=expression("E(N|"~(theta~"',"~ theta~"'")~")" ) 
name3=expression("E(N|"~(theta~"',"~ 0)~")" ) 
name4=expression("E(N|"~(theta~"',"~ infinity)~")" ) 
name5=expression("E(N|"~(0~","~ theta~"'")~")" ) 
name6=expression("E(N|"~(0~","~ 0)~")" ) 
name7=expression("E(N|"~(infinity ~","~ theta~"'")~")" ) 


legend(130, 465, legend=c(name1, name2, name3, name4,name5,name6,name7 ),
       col=c("Black","Dark Blue","Dark Green","Purple","Orange","Red","Blue"),lwd = 2, cex=0.8,
       box.lty=0)

points(x=0:196,y=ESSeachrunTPall200$line1[1:197],type="l",col="Dark Blue",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line2[1:197],type="l",col="Dark Green",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line3[1:197],type="l",col="Purple",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line4[1:197],type="l",col="Orange",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line5[1:197],type="l",col="Red",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line6[1:197],type="l",col="Blue",lwd = 2,lty=1)


name2=expression(alpha[sub1^"'"*sub2])

load(file="SamplesizesepateTP.RData")

SamplesizesepateTP$MaxsepateTP
#abline(h=SamplesizesepatePWP$MaxsepatePWP,lty=3,col="black")
#abline(h=SamplesizesepatePWP$Line1sepatePWP,lty=3,col="Dark Blue")
#need to find when this goes above


nopointmax=min(which(SamplesizesepateTP$MaxsepateTP<=maxeachrunTP200))-1
nopointline1=min(which(SamplesizesepateTP$Line1sepateTP<=ESSeachrunTPall200$line1))-1
nopointline2=min(which(SamplesizesepateTP$Line2sepateTP<=ESSeachrunTPall200$line2))-1
nopointline3=min(which(SamplesizesepateTP$Line3sepateTP<=ESSeachrunTPall200$line3))-1
nopointline4=min(which(SamplesizesepateTP$Line4sepateTP<=ESSeachrunTPall200$line4))-1
nopointline5=min(which(SamplesizesepateTP$Line5sepateTP<=ESSeachrunTPall200$line5))-1
nopointline6=min(which(SamplesizesepateTP$Line6sepateTP<=ESSeachrunTPall200$line6))-1


abline(v=nopointmax,lty=3,lwd = 3,col="Black")
abline(v=nopointline1,lty=3,lwd = 3,col="Dark Blue")
abline(v=nopointline2,lty=3,lwd = 3,col="Dark Green")
abline(v=nopointline3,lty=3,lwd = 3,col="Purple")
abline(v=nopointline4,lty=3,lwd = 3,col="Orange")
abline(v=nopointline5,lty=3,lwd = 3,col="Red")
abline(v=nopointline6,lty=3,lwd = 3,col="Blue")
#4.62*1.2 by 5.6*1.2 
#5.54 by 6.72
#TPplotSSR2
#5.54 by 6.72
#TPplotS2

#values for paper 
min(maxeachrunTP200[1:197])
max(maxeachrunTP200[1:197])
nopointmax
SamplesizesepateTP$MaxsepateTP

SamplesizesepateTP$Line6sepateTP
nopointline6

#FWER control version of all the plots


#pairwisepower
load(file = "ESSeachrunPWPall200.RData")
load(file="maxeachrunPWP200.RData")
minPW=min(ESSeachrunPWPall200$line5[1:155]) #minimum value
maxPW=max(maxeachrunPWP200[1:155]) #minimum value
length(maxeachrunPWP200[1:155])
#plot
plot(x=0:154,y=maxeachrunPWP200[1:155],xlim=c(0,154),ylim = c(minPW-50,maxPW+10) ,main="Setting 1",type="l",col="Black",lty=1,lwd = 2,xlab = "When the second active treatment is added",ylab = "Required sample size")

#name1=expression("E(N|" (theta[1] == 0)) expression((theta[1] == 0))

name1=expression("max(N)")
name2=expression("E(N|"~(theta~"',"~ theta~"'")~")" ) 
name3=expression("E(N|"~(theta~"',"~ 0)~")" ) 
name4=expression("E(N|"~(theta~"',"~ -infinity)~")" ) 
name5=expression("E(N|"~(0~","~ theta~"'")~")" ) 
name6=expression("E(N|"~(0~","~ 0)~")" ) 
name7=expression("E(N|"~(-infinity ~","~ theta~"'")~")" ) 


legend(100, 370, legend=c(name1, name2, name3, name4,name5,name6,name7 ),
       col=c("Black","Dark Blue","Dark Green","Purple","Orange","Red","Blue"),lwd = 2, cex=0.8,
       box.lty=0)

points(x=0:154,y=ESSeachrunPWPall200$line1[1:155],type="l",col="Dark Blue",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line2[1:155],type="l",col="Dark Green",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line3[1:155],type="l",col="Purple",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line4[1:155],type="l",col="Orange",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line5[1:155],type="l",col="Red",lwd = 2,lty=1)
points(x=0:154,y=ESSeachrunPWPall200$line6[1:155],type="l",col="Blue",lwd = 2,lty=1)


name2=expression(alpha[sub1^"'"*sub2])

load(file="SamplesizesepatePWPFWER.RData")

SamplesizesepatePWPFWER$MaxsepatePWP
#abline(h=SamplesizesepatePWP$MaxsepatePWP,lty=3,col="black")
#abline(h=SamplesizesepatePWP$Line1sepatePWP,lty=3,col="Dark Blue")
#need to find when this goes above


nopointmax=min(which(SamplesizesepatePWPFWER$MaxsepatePWP<=maxeachrunPWP200))-1
nopointline1=min(which(SamplesizesepatePWPFWER$Line1sepatePWP<=ESSeachrunPWPall200$line1))-1
nopointline2=min(which(SamplesizesepatePWPFWER$Line2sepatePWP<=ESSeachrunPWPall200$line2))-1
nopointline3=min(which(SamplesizesepatePWPFWER$Line3sepatePWP<=ESSeachrunPWPall200$line3))-1
nopointline4=min(which(SamplesizesepatePWPFWER$Line4sepatePWP<=ESSeachrunPWPall200$line4))-1
nopointline5=min(which(SamplesizesepatePWPFWER$Line5sepatePWP<=ESSeachrunPWPall200$line5))-1
nopointline6=min(which(SamplesizesepatePWPFWER$Line6sepatePWP<=ESSeachrunPWPall200$line6))-1


abline(v=nopointmax,lty=3,lwd = 3,col="Black")
abline(v=nopointline1,lty=3,lwd = 3,col="Dark Blue")
abline(v=nopointline2,lty=3,lwd = 3,col="Dark Green")
abline(v=nopointline3,lty=3,lwd = 3,col="Purple")
abline(v=nopointline4,lty=3,lwd = 3,col="Orange")
abline(v=nopointline5,lty=3,lwd = 3,col="Red")
abline(v=nopointline6,lty=3,lwd = 3,col="Blue")
#5.54 by 6.72
#PWPplotSSR2
#PWPplotS1
#values for paper 
min(maxeachrunPWP200)
max(maxeachrunPWP200)
nopointmax
SamplesizesepatePWPFWER$MaxsepatePWP

nopointline6
SamplesizesepatePWPFWER$Line6sepatePWP

#totalpower
load(file = "ESSeachrunTPall200.RData")
load(file="maxeachrunTP200.RData")
minTP=min(ESSeachrunTPall200$line5[1:197]) #minimum value
maxTP=max(maxeachrunTP200[1:197]) #minimum value
length(maxeachrunTP200[1:197])
#plot
plot(x=0:196,y=maxeachrunTP200[1:197],xlim=c(0,196),ylim = c(minTP-50,maxTP+10) ,main="Setting 1" ,type="l",col="Black",lty=1,lwd = 2,xlab = "When the second active treatment is added",ylab = "Required sample size")

#name1=expression("E(N|" (theta[1] == 0)) expression((theta[1] == 0))

name1=expression("max(N)")
name2=expression("E(N|"~(theta~"',"~ theta~"'")~")" ) 
name3=expression("E(N|"~(theta~"',"~ 0)~")" ) 
name4=expression("E(N|"~(theta~"',"~ -infinity)~")" ) 
name5=expression("E(N|"~(0~","~ theta~"'")~")" ) 
name6=expression("E(N|"~(0~","~ 0)~")" ) 
name7=expression("E(N|"~(-infinity ~","~ theta~"'")~")" ) 


legend(130, 465, legend=c(name1, name2, name3, name4,name5,name6,name7 ),
       col=c("Black","Dark Blue","Dark Green","Purple","Orange","Red","Blue"),lwd = 2, cex=0.8,
       box.lty=0)

points(x=0:196,y=ESSeachrunTPall200$line1[1:197],type="l",col="Dark Blue",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line2[1:197],type="l",col="Dark Green",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line3[1:197],type="l",col="Purple",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line4[1:197],type="l",col="Orange",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line5[1:197],type="l",col="Red",lwd = 2,lty=1)
points(x=0:196,y=ESSeachrunTPall200$line6[1:197],type="l",col="Blue",lwd = 2,lty=1)


name2=expression(alpha[sub1^"'"*sub2])

load(file="SamplesizesepateTPFWER.RData")

SamplesizesepateTPFWER$MaxsepateTP



nopointmax=min(which(SamplesizesepateTPFWER$MaxsepateTP<=maxeachrunTP200))-1
nopointline1=min(which(SamplesizesepateTPFWER$Line1sepateTP<=ESSeachrunTPall200$line1))-1
nopointline2=min(which(SamplesizesepateTPFWER$Line2sepateTP<=ESSeachrunTPall200$line2))-1
nopointline3=min(which(SamplesizesepateTPFWER$Line3sepateTP<=ESSeachrunTPall200$line3))-1
nopointline4=min(which(SamplesizesepateTPFWER$Line4sepateTP<=ESSeachrunTPall200$line4))-1
nopointline5=min(which(SamplesizesepateTPFWER$Line5sepateTP<=ESSeachrunTPall200$line5))-1
nopointline6=min(which(SamplesizesepateTPFWER$Line6sepateTP<=ESSeachrunTPall200$line6))-1


abline(v=nopointmax,lty=3,lwd = 3,col="Black")
abline(v=nopointline1,lty=3,lwd = 3,col="Dark Blue")
abline(v=nopointline2,lty=3,lwd = 3,col="Dark Green")
abline(v=nopointline3,lty=3,lwd = 3,col="Purple")
abline(v=nopointline4,lty=3,lwd = 3,col="Orange")
abline(v=nopointline5,lty=3,lwd = 3,col="Red")
abline(v=nopointline6,lty=3,lwd = 3,col="Blue")

