setwd("E:/hjl/bubble")
testall <- read.csv("testall2.csv", header = F)
testall[, 1] <- as.Date(testall[, 1])
par(mar = c(2.4,2.5,1,0)+0.1)
plot(testall, type = "l", lwd = 1.5, ylab = "", xlab = "", cex.axis = .8 )
abline(v = testall[c(6, 18, 62, 72, 108,126,173, 203, 215, 225), 1], lty = 2)
polygon(testall[c(6:18, 18:6), 1], c(rep(7000, length(6:18)), 
rep(-500, length(6:18))),col = "grey", lty =2)
polygon(testall[c(62:72, 72:62),1], c(rep(7000, length(62:72)), 
rep(-500, length(62:72))), col = "grey", lty =2)
polygon(testall[c(108:126, 126:108),1], c(rep(7000, length(108:126)), 
rep(-500, length(108:126))), col = "grey", lty =2)
polygon(testall[c(173:203, 203:173),1], c(rep(7000, length(173:203)), 
rep(-500, length(173:203))), col = "grey", lty =2)
polygon(testall[c(215:225, 225:215),1], c(rep(7000, length(215:225)), 
rep(-500, length(215:225))), col = "grey", lty =2)
lines(testall, type = "l", lwd = 1.5, ylab = "", xlab = "")


text(testall[11,1], 2000, "ADF = 6.252", cex= .8)
text(testall[68,1], 3000, "ADF = 3.19", cex= .8)
text(testall[117,1], 4000, "ADF = 1.83", cex= .8)
text(testall[190,1], 5000, "ADF = 3.58", cex= .8)
text(testall[220,1], 5000, "ADF = 4.129", cex= .8)
text(testall[2, 1], 6000, "1991/05/31", cex = .8)
text(testall[22, 1], 6000, "1992/05/29", cex = .8)
text(testall[57, 1], 6000, "1996/01/31", cex = .8)
text(testall[77, 1], 6000, "1996/11/29", cex = .8)
text(testall[106, 1], 6000, "1999/11/30", cex = .8)
text(testall[128, 1], 6000, "2001/05/31", cex = .8)
text(testall[171, 1], 6000, "2005/04/29", cex = .8)
text(testall[194, 1], 6000, "2007/10/31", cex = .8)
text(testall[213, 1], 6000, "2008/10/31", cex = .8)
text(testall[233, 1], 6000, "2009/08/31", cex = .8)

##
test <- function(T) {
	mu <- 0.0024
	sigmad <- 0.001
	D_0 <- 1
	rho <- 0.985
	b <- 1
	B <- vector()
	B[1] <- 0.5
	pi <- 0.85
	eta <- 0.5
	tao <- 0.05
	k <- 50
	y_t <- rnorm(T, sd = tao)
	e_bt <- exp(y_t - tao^2/2)
	theta <- rbinom(T, 1, pi)
	for (i in 2:T){
		if(B[i - 1] < b)
		B[i] <- rho^(-1) * B[i - 1] * e_bt[i]
		else
		B[i] <- (eta + (pi * rho)^(-1) * theta[i] * (B[i - 1] - rho * eta)) * e_bt[i]
	}
	# return(B = B + cumsum(rnorm(T)))

	rho0 <- 0.9
	temp <- vector()
	temp[1] <- rnorm(1)
	for(ii in 2:T)
		temp[ii] < rho0 * temp[ii - 1] + rnorm(1)
	return(B = B + temp)	 
}



library(fUnitRoots)

powerSim <- function(T, n) {
adfunit <- rep(0, n)
pwy <- rep(0, n)
hg <- rep(0, n)

for (i in 1:n) {	
	dat <- test(T)
	# plot(dat, type = "l")
	adfunit[i] <- ifelse(adfTest(dat)@test$stat > unidf[T], 1, 0)
	for(j in 20:T) {
		if(adfTest(dat[1:j])@test$stat > pwydf[T])
		{pwy[i] <- 1
		break}			
	}

	for (k in 1:(T - 20)) {
		for(m in (k + 20):T) {
		if(adfTest(dat[k:m])@test$stat > hgdf[T])
		{hg[i] <- 1
		break}
		}
		if(hg[i] == 1)
		break
	}	
	cat(i,'\n');
}
return(c(mean(adfunit),mean(pwy),mean(hg)));
}

##
test2 <- function(T) {
	mu <- 0.0024
	sigmad <- 0.001
	D_0 <- 1
	rho <- 0.985
	b <- 1
	B <- vector()
	B[1] <- 0.5
	pi <- 0.85
	eta <- 0.5
	tao <- 0.05
	k <- 50
	flag=0;
	while (flag<2)
	{
	bk=0;
	y_t <- rnorm(T, sd = tao)
	e_bt <- exp(y_t - tao^2/2)
	theta <- rbinom(T, 1, pi)
	for (ii in 2:T){
	  if(bk<2)
		if(B[ii - 1] < b)
		B[ii] <- rho^(-1) * B[ii - 1] * e_bt[ii]
		else
		{
		B[ii] <- (eta + (pi * rho)^(-1) * theta[ii] * (B[ii - 1] - rho * eta)) * e_bt[ii]
		if (theta[ii]==0) bk=bk+1;
		}
	  else
        B[ii]=B[ii-1]*e_bt[ii];  	  
	}
	flag=bk;
	}
	# return(B = B + cumsum(rnorm(T)))

	rho0 <- 0.9
	temp <- vector()
	temp[1] <- rnorm(1)
	for(ii in 2:T)
		temp[ii] < rho0 * temp[ii - 1] + rnorm(1)
	return(B = B + temp)	 
}

##test for two bubbles
library(fUnitRoots)

powerSim2 <- function(T, n) {
pwy1 <- rep(0, n)
hg1 <- rep(0, n)
pwy <- rep(0, n)
hg <- rep(0, n)
for (i in 1:n) {	
	dat <- test2(T)
	# plot(dat, type = "l")
	f1=0;f2=0;rem=0;
	hga=matrix(-100,T,T);
	for(j in 20:T) {
       tstat=adfTest(dat[1:j])@test$stat;	
	   if (f1==0){
		if(tstat > pwydf[T])
		{f1= 1;
		rem=j;
		pwy1[i]=1;}
		}
	   if(f1==1)
	   {
		if (f2==0){
		if ((tstat<pwydf[T])&(j>rem+10)) f2=1
		}
		if ((f2==1)&(tstat > pwydf[T])) {
		pwy[i]=1;
		break
        }	
       }		
	}
    maxx=-100; 
	for (k in 1:(T - 20)) {
		for(m in (k + 20):T) {
		hga[k,m]=adfTest(dat[k:m])@test$stat
		if (hga[k,m]>maxx) {maxx=hga[k,m];maxi=k;maxj=m;}
		}
	}
	if (maxx>hgdf[T])
	{
	hg1[i]=1;
	if (maxi==1) tempi=1 else tempi=maxi-1;
	if (maxj==T) tempj=T else tempj=maxj+1;
 	if (max(max(hga[1:tempi,1:tempi]),max(hga[tempj:T,tempj:T]))>hgdf[T]) hg[i]=1;
	}
	cat(i,'\n');
}
return(c(mean(pwy1),mean(hg1),mean(pwy),mean(hg)));
}

T = c(100, 200, 400)
n = 500
pwydf=rep(0,400)
pwydf[100]=1.078;
pwydf[200]=1.229;
pwydf[400]=1.341;
hgdf=rep(0,400);
hgdf[100]=1.489;
hgdf[200]=1.730;
hgdf[400]=1.893;
unidf=rep(0,400);
unidf[100]=-0.09;
unidf[200]=-0.061;
unidf[400]=-0.075;
a1=powerSim2(T[1], n)
a2=powerSim2(T[2], n)
a3=powerSim2(T[3], n)
write.csv(c(a1,a2,a3),"result2.csv");
b1=powerSim(T[1], n)
b2=powerSim(T[2], n)
b3=powerSim(T[3], n)
write.csv(c(b1,b2,b3),"onebubblepower1.csv");
library(xtable)

data2=test2(400);
plot(data2,type='l');
