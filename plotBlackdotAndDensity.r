
## iteration 40000 State 2452 alpha 9183 rho 16006 Bsigma 6889 sigma 1748 iterkappa 2  

postate<-read.csv(file="KappaState2019-05-19.txt",sep="",header=FALSE)

end=max(dim(postate) )
burnin = 10000 #30000
mstate<-apply(postate[burnin:end,],2,mean)



pdf('KappaStateup1point05-19.pdf')
plot(ss[, 1],col="red",main="animal 1",xlab="times index",ylab="state",ylim=c(1,2) )
points(mstate[1:50],pch="+")
abline(h=1.5)
plot(ss[, 2],col="red",main="animal 2",xlab="times index",ylab="state",ylim=c(1,2) )
points(mstate[51:100],pch="+")
abline(h=1.5)
plot(ss[, 3],col="red",main="animal 3",xlab="times index",ylab="state",ylim=c(1,2) )
points(mstate[101:150],pch="+")
abline(h=1.5)
plot(ss[, 4],col="red",main="animal 4",xlab="times index",ylab="state",ylim=c(1,2) )
points(mstate[151:200],pch="+")
abline(h=1.5)
plot(ss[, 5],col="red",main="animal 5",xlab="times index",ylab="state",ylim=c(1,2) )
points(mstate[201:250],pch="+")
abline(h=1.5)
dev.off()




llambda<-read.csv(file="KappaLambda2019-05-19.txt",sep="",header=FALSE)
mean(llambda[burnin:end,1])
sd(llambda[burnin:end,1])

lalpha<-read.csv(file="blacBGalpha2019-05-19.txt",sep="",header=FALSE)
#burnin = 5000; end = max( dim(lalpha) )
mean(lalpha[burnin:end,1])
sd(lalpha[burnin:end,1])


lrho<-read.csv(file="blacBGrho2019-05-19.txt",sep="",header=FALSE)
mean(lrho[burnin:end,1])
sd(lrho[burnin:end,1])

lBsig<-read.csv(file="blacBGBsig2019-05-19.txt",sep="",header=FALSE)
mean(lBsig[burnin:end,1])
sd(lBsig[burnin:end,1])

lsigma<-read.csv(file="blacBGsigma2019-05-19.txt",sep="",header=FALSE)
mean(lsigma[burnin:end,1])
sd(lsigma[burnin:end,1])


alpha <- 1.2 #0.8
rho   <- 5  #5.2
sigma <- 0.7   # 2.7
Bsigma <- 2
swlamda = c( 0.1, 0.4)


# bw.nrd(lalpha[burnin:end,1])  ## choose bandwidth  bw 
dal=density(lalpha[burnin:end,1],bw=0.02,kernel="gaussian")
plot(dal,type="n",main=expression("(c) Posterior for"~alpha),xlab="")
polygon(dal, col="#ff606080", border=NA)
abline(v=alpha,col='blue',lwd=2,lty=2)


drho = density(lrho[burnin:end,1],bw=0.1,kernel="gaussian")
plot(drho,type='n',main=expression("(d) Posterior for"~rho),xlab="")
polygon(drho, col="#ff606080", border=NA)
abline(v=rho,col='blue',lwd=2,lty=2)


dsig=density(lsigma[burnin:end,1],bw=0.02,kernel="gaussian")
plot(dsig,type='n',main=expression("(f) Posterior for"~sigma),xlab="")
polygon(dsig, col="#ff606080", border=NA)
abline(v=sigma,col='blue',lwd=2,lty=2)


# bw.nrd(lBsig[burnin:end,1])
dBsig=density(lBsig[burnin:end,1],bw=0.05,kernel="gaussian")
plot(dBsig,type='n',main=expression("(f) Posterior for"~sigma),xlab="")
polygon(dBsig, col="#ff606080", border=NA)
abline(v=Bsigma,col='blue',lwd=2,lty=2)



ldlam1=density(llambda[burnin:end,1],bw=0.04,kernel="gaussian")
#usex1 = ldlam1$x[which(ldlam1$x>0&ldlam1$x<0.02)]
#usey1 = ldlam1$y[which(ldlam1$x>0&ldlam1$x<0.02)]
usex2 = ldlam1$x[which(ldlam1$x>0.3)]
usey2 = ldlam1$y[which(ldlam1$x>0.3)]

xlam = llambda[burnin:end,1]
h <- density(xlam, kernel="gaussian")$bw
w <- 1 / pnorm(0, mean=xlam, sd=h, lower.tail=FALSE)
# bw.nrd(llambda[burnin:end,1])
dlam=density(xlam,bw=0.04,kernel="gaussian", weights=w / length(xlam),cut=0)
# dlam=density(llambda[burnin:end,1],bw=0.04,kernel="gaussian")
# dlam$x = c(0,usex1,dlam$x,usex2)  ;  dlam$y = c(0,usey1,dlam$y,usey2)
dlam$x = c(0,dlam$x,usex2)  ;  dlam$y = c(0,dlam$y,usey2)
plot(dlam,type='n',main=expression(paste("(f) Posterior for ",lambda,"1")),xlab="")#,xlim=c(-0.01,0.4))#,ylim=c(0.03,0.2) )
polygon(c(dlam$x),c(dlam$y), col="#ff606080", border=NA)
#polygon(dlam, col="#ff606080", border=NA)
abline(v=swlamda[1],col='blue',lwd=2,lty=2)



ldlam2=density(llambda[burnin:end,2],bw=0.08,kernel="gaussian")
usex = ldlam2$x[which(ldlam2$x>0.6)]
usey = ldlam2$y[which(ldlam2$x>0.6)]

xlam2 = llambda[burnin:end,2]
h <- density(xlam2, kernel="gaussian")$bw
w <- 1 / pnorm(0, mean=xlam2, sd=h, lower.tail=FALSE)
dlam2=density(xlam2,bw=0.08,kernel="gaussian", weights=w / length(xlam2),cut=0)
dlam2$x = c(0,dlam2$x,usex)  ;  dlam2$y = c(0,dlam2$y,usey)

plot(dlam2,type='n',main=expression(paste("(f) Posterior for ",lambda,"2")),xlab="")
polygon(c(dlam2$x),c(dlam2$y), col="#ff606080", border=NA)
#polygon(dlam2, col="#ff606080", border=NA)
abline(v=swlamda[2],col='blue',lwd=2,lty=2)




pdf('density190519.pdf')
par(mfrow=c(3,2))
plot(dal,type="n",main=expression("(c) Posterior for"~alpha),xlab="")
polygon(dal, col="#ff606080", border=NA)
abline(v=alpha,col='blue',lwd=2,lty=2)

plot(drho,type='n',main=expression("(d) Posterior for"~rho),xlab="")
polygon(drho, col="#ff606080", border=NA)
abline(v=rho,col='blue',lwd=2,lty=2)

plot(dsig,type='n',main=expression("(f) Posterior for"~sigma),xlab="")
polygon(dsig, col="#ff606080", border=NA)
abline(v=sigma,col='blue',lwd=2,lty=2)

plot(dBsig,type='n',main=expression("(f) Posterior for BM"~sigma),xlab="")
polygon(dBsig, col="#ff606080", border=NA)
abline(v=Bsigma,col='blue',lwd=2,lty=2)

plot(dlam,type='n',main=expression(paste("(f) Posterior for",lambda,"1")),xlab="")#,xlim=c(-0.01,0.4))#,ylim=c(0.03,0.2) )
polygon(c(dlam$x),c(dlam$y), col="#ff606080", border=NA)
#polygon(dlam, col="#ff606080", border=NA)
abline(v=swlamda[1],col='blue',lwd=2,lty=2)

plot(dlam2,type='n',main=expression(paste("(f) Posterior for ",lambda,"2")),xlab="")
polygon(c(dlam2$x),c(dlam2$y), col="#ff606080", border=NA)
abline(v=swlamda[2],col='blue',lwd=2,lty=2)

dev.off()








