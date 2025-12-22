
#set working directory and source files
setwd("D:/1. Yen/1. Study/2. Thesis")
source("sviRoots.R")					#includes function for crossedness by Jim Gatheral (the code and his lectures are available online)
source("Functions built for calibration.R")				#includes all self-defined functions

#1. IMPORT AND HANDLE DAT-FILE DATA: from alloptions->newdata
#2. F AND DF <-PUT CALL PARITY
#3. SEE THE PLOT IV SMILE
#4. INITIAL GUESS
#5. OPTIMAL 

#1.IMPORT AND HANDLE DAT-FILE DATA: from rawdata->newdata
#1.1 import file to R

#  import dataset name"quotedata". then drop the first 2 rows
rawdata=quotedata[c(-1:-3),c(-9,-10,-20,-21)] #drop delta, gamma
colnames(rawdata)=c("Expirationdate","call","lastsale", "net","bid", "ask", "vol", "iv", "openint", "strike","put","lastsale", "net","bid", "ask", "vol", "iv", "openint");
rawdata$call <- "C" # change call options to "C"
rawdata$put <-"P"
colnames(rawdata)=c("Expirationdate","optiontype","lastsale", "net","bid", "ask", "vol", "iv", "openint", "strike","optiontype","lastsale", "net","bid", "ask", "vol", "iv", "openint");

# Expiration Date,Calls,Last Sale,Net,Bid,Ask,Vol(volum) ,IV,Open Int,Strike,
# Puts,Last Sale,Net,Bid,Ask,Vol(volum) ,IV,Open Int.

#format data
spot= as.numeric(quotedata[1,2])
quotedate <- as.Date("05/13/19",format = "%m/%d/%y") #format the start date 
rawdata$Expirationdate=as.Date(rawdata$Expirationdate, format = "%m/%d/%Y") #format the end date 

#add columns: quote date, maturity
rawdata$quotedate=quotedate
rawdata$ttm=round(as.numeric((rawdata$Expirationdate - rawdata$quotedate)/365), digits = 2)

calloptions=rawdata[ ,c(19,1,20,10,2:9)] #re-arrange the order of the columns
putoptions=rawdata[ ,c(19,1,20,10,11:18)]
alloptions=rbind(calloptions,putoptions) #combine call then put into a dataframe

alloptions$strike= as.numeric(alloptions$strike)
alloptions$bid=as.numeric(alloptions$bid)
alloptions$ask=as.numeric(alloptions$ask)
alloptions$iv=as.numeric(alloptions$iv)
alloptions$vol=as.numeric(alloptions$vol) #vol is volume

alloptions$pricemid=0.5*(alloptions$bid+alloptions$ask)

#1.2 Data handle

#choose: bid, ask>0; volume> 0, IV >0, maturity: 5 days - 2 years
newdata <- alloptions[which(alloptions$bid>0 & alloptions$ask > 0 & alloptions$vol>0 & 
                              alloptions$iv>0 & alloptions$ttm > 5/365 & alloptions$ttm < 2), ]
allmaturities=unique(newdata$ttm) #see all maturities. we see that we have 28 maturities.
selected.maturities = c(0.03, 0.07, 0.11, 0.18, 0.36, 0.61, 1.1, 1.6) #choose 8 out of 28 maturities
numslices=8
selected.data=newdata[which(newdata$ttm %in% selected.maturities), ]
#selected.data = selected.data[ , c(-11)] #drop iv

#then create 3 data sets: 1. data.complete, 2. DF_Flist, 3. dataset.unique

data.complete <- data.frame(matrix(ncol = 15, nrow = 0)) #create an empty data complete. (here strike is duplicated)
colnames(data.complete ) <- c("quotedate","Expirationdate","ttm","strike", "optiontype","lastsale", "net","bid", "ask", 
                              "vol", "iv", "openint","pricemid", "DF", "forward")
DF_Flist <- data.frame(matrix(ncol = 3, nrow = 0)) #create a list of DF and forward
colnames(DF_Flist) <- c("ttm","DF","forward")

dataset.unique= data.frame(matrix(ncol = 7, nrow = 0)) #create data.unique. CALL DATA. (here strike is NOT duplicated)
colnames(dataset.unique) <- c("ttm","strike","iv","pricemid", "DF","forward", "ivol")
#ivol is implied vol from (mkt opt proce and BS formula)
#iv in quotedata: what is it?

#2. DF and forward
for (year in selected.maturities){ #One expiration at a time.

  vdc <- selected.data[(selected.data$ttm==year)& selected.data$optiontype=='C',];# Select just this expiration
  vdp <- selected.data[(selected.data$ttm==year)& selected.data$optiontype=='P',];# Select just this expiration
  callStrikes <- sort(unique(vdc$strike));
  putStrikes <- sort(unique(vdp$strike));
  strikes <- callStrikes[callStrikes%in%putStrikes];
 
  decoratedoptions=selected.data[(selected.data$ttm==year & selected.data$strike%in%strikes), ]
  #one K may occur twice
  
  nK <- length(strikes); #the 2 section below are for C-P parity
  ca <- numeric(nK);
  cb <- numeric(nK);
  pa <- numeric(nK);
  pb <- numeric(nK);
  imid <- numeric(nK);
  
  cbb <- vdc$bid; #call bid bid 
  pbb <- vdp$bid;
  cba <- vdc$ask;
  pba <- vdp$ask;
  
  civ <- vdc$iv; #iv. #this section is for taking call data. unique (no duplicated K)
  cpricemid<- vdc$pricemid
  iv <- numeric(nK) #implied vol
  pricemid <-numeric(nK) #price mid
    
    for (i in 1:nK){
      K <- strikes[i]; #for each K, we might have many call bid prices --> call them: cal bid bid (cbb)
      
      cb[i] <- mean(cbb[vdc$strike==K]); #call bid = mean(call bid bid)
      pb[i] <- mean(pbb[vdp$strike==K]); #put bid 
      ca[i] <- mean(cba[vdc$strike==K]); #call ask
      pa[i] <- mean(pba[vdp$strike==K]); #put ask
      
      ibid <- cb[i]-pb[i]; # (c-P) bid
      iask <- ca[i]-pa[i];
      imid[i] <- (ibid+iask)/2; #(C-P) mid
      
      iv[i] <- mean(civ[vdc$strike==K]); #iv= mean(iv)
      pricemid[i]<- mean(cpricemid[vdc$strike==K])
      }
 
  x=strikes
  y=imid
  lm=lm(y~x,data.frame(x,y))   #Create the linear regression
  plot(y~x,data = data.frame(x,y), xlab = "strike",ylab = "CPdifference", main = "linear regression") #Plot the results
  abline(lm) #Add a regression line
  summary(lm)   #view the result. pvalue should be smaller than alpha=0.1
  
  pvalues=summary(lm)$coefficients[,4] 
  pvalue1= pvalues["(Intercept)"]
  pvalue2= pvalues["x"]
  
  if (is.na(pvalue1) || is.na(pvalue2)){ #there is no pair of CP --> no regression
    print("linear regression failed") 
  } else if (pvalue1 >0.01 || pvalue2 > 0.01) {
    print("the result of linear regression is not significant") 
  
  } else if  (pvalue1 <0.01 & pvalue2 < 0.01){ # pvalue should be smaller than alpha=0.1
    intercept= as.numeric(coef(lm)["(Intercept)"])
    slope = as.numeric(coef(lm)["x"])
    DF = -slope
    forward = intercept/DF
    linearFormula = intercept+slope*strikes #linearformula = expected C-P (from the regression result)    
    error = (imid-linearFormula)/forward #(actual y - expected y)/forward: normalize things
    errorSlice = mean(error) 
    
    if (errorSlice > 0.01){
      print("DF and forward extraction failed")
    } else {
    decoratedoptions$DF = DF    #this is for data complete
    decoratedoptions$forward = forward
    
    DF_FoneT=data.frame("ttm"=year,"DF"=DF, "forward"=forward) #this is for DF_F list
    
    datasetoneT.unique=data.frame("ttm"=year,"strike"=strikes, "iv"=iv, "pricemid"=pricemid,"DF"=DF, "forward"=forward)
    datasetoneT.unique$ivol= BSImpliedVolCall(year, strikes, DF, forward, pricemid)
   #check: BSImpliedVolCall(year, strike=2890, DF=0.9985929, forward=2850.018, 22.700) calulated: 0.1544864, table: 0.1544869
    } 
  
  } #end of the 1st if

  data.complete=rbind(data.complete,decoratedoptions) #this is for data complete. all data for 6 selected maturities
  DF_Flist = rbind(DF_Flist, DF_FoneT) ##this is for DF_F list. list of DF and forward
  dataset.unique=rbind(dataset.unique, datasetoneT.unique) #this is for taking call data. unique (no duplicated K)
  
}#end of for
dataset= data.complete[data.complete$optiontype=="C", ]#final dataset: all CALL data for 6T, but duplicated K
colnames(dataset.unique)=c("ttm", "strike","iv","pricemid", "DF", "forward", "ivol")


#3. SEE THE PLOT IV SMILE
plot(log(datasetoneT.unique$strike/datasetoneT.unique$forward),datasetoneT.unique$ivol, xlab ="Log-Strike", 
     ylab ='Implied volatility', main= 'T=1.6') #iv smile


#4. INITIAL GUESS 
#4.1 ATM Total IV computation: use dataset.unique
#4.2 inital guess
#4.3 plot the fit

#4.1 ATM Total IV computation: use dataset.unique --> w0
library("stinepack")
w0 = numeric(numslices)
for (i in 1:numslices){
  
  t = DF_Flist$ttm[i]
  impliedvol = dataset.unique[which(dataset.unique$ttm==t), ]$ivol
  midvar = (impliedvol)^2
  forward = dataset.unique[which(dataset.unique$ttm==t), ]$forward
  strike = dataset.unique[which(dataset.unique$ttm==t), ]$strike
  k = log(strike/forward)
  w0 [i] = t*(stinterp(k,midvar,0)$y)
  if (is.na(w0[i])==TRUE){
  w0 [i] = t* (midvar[which.min(abs(k))])  
  }
}
DF_Flist$w0= w0

#combine w0 with dataset.unique
dataset.unique$w0= 0
for (year in dataset.unique$ttm){
  dataset.unique$w0 [which(dataset.unique$ttm==year)]=DF_Flist$w0 [which(DF_Flist$ttm==year)]
}
dataset.unique=dataset.unique[c(1,2,5,6,8,4,3,7)]

#4.2 initial guess
plot(DF_Flist$ttm,DF_Flist$w0) #see plot of w0
library("GenSA")
#simulated annealing
out <- GenSA(par=c(eta = 1, rho = 0.5), fn=SSVIpriceSA , lower = c(eta=0, rho=-1), upper = c(eta=10, rho=1), 
             control=list(maxit=500, smooth=FALSE), datavec=dataset.unique); 
                                #maxit: Integer. Max number of iterations of the algorithm
out[c("value","par","counts")]; #value: value of fn corresponding to par. 
                                #par: Vector. The best set of parameters found. 
out$par[1]*(1+abs(out$par[2]));		#check if constraint eta*(1+|rho|)<2 is met. the result MUST < 2

#fitting total implied variance surface versus data
params <- out$par; #then remove out
surfacesliceparams <- array(NA, dim = c(numslices, 5));
atmtvs <- DF_Flist$w0	
ttm= DF_Flist$ttm
for (i in 1:numslices){
  paramsJW <- SSVItoJW(params[1],params[2], atmtvs[i] ,ttm[i]);
  surfacesliceparams[i,] <- JWtoRAW(c(paramsJW,ttm[i]));
}
colnames(surfacesliceparams)= c("a","b",'rho','m', "sig")

#4.3 PLOT 2D OF THE FIT
summederror <- 0
summedcrossedness <- 0 
slicecol <- rainbow(numslices)
for (i in 1:numslices){
 
  t = DF_Flist$ttm[i]
  dataset_slice_SA = dataset.unique[which(dataset.unique$ttm==t), ];
  slice_k = log(dataset_slice_SA$strike/dataset_slice_SA$forward);	
  slice_w = ((dataset_slice_SA$ivol)^2)*t
  summederror <- summederror + wrawpriceSA(surfacesliceparams[i,], dataset_slice_SA);
  if(i==1){
    plot(slice_k, slice_w, col = slicecol[i],type="p", xlab = "Log moneyness k",  xlim=c(-0.45,0.2), ylim = c(0,0.08),
         ylab = "Total implied variance w(k,T)", main="SSVI fit to volatility surface" )
  } else{
    points(slice_k, slice_w, col = slicecol[i], type="p");
    crossedness_arg = as.data.frame(rbind(surfacesliceparams[i-1,c(1,2,5,3,4)], surfacesliceparams[i,c(1,2,5,3,4)]));
    colnames(crossedness_arg)<- c("a","b","sig","rho","m");
    summedcrossedness <- summedcrossedness + ComputeSviRoots(crossedness_arg)$crossedness;
  }
  slicesize <- 100;
  slicefitdata <- array(NA, dim=c(slicesize,7));							
  slicefitdata[,1] <- seq(from = -2, to = 1, length.out = slicesize);			#log moneyness k
  slicefitdata[,2] <- rep(surfacesliceparams[i,1], times = slicesize);		#a
  slicefitdata[,3] <- rep(surfacesliceparams[i,2], times = slicesize);    #b
  slicefitdata[,4] <- rep(surfacesliceparams[i,3], times = slicesize);    #rho
  slicefitdata[,5] <- rep(surfacesliceparams[i,4], times = slicesize);    #m
  slicefitdata[,6] <- rep(surfacesliceparams[i,5], times = slicesize);    #sigma
  slicefitdata[,7] <- apply(slicefitdata[,1:6], 1, wraw2);                #w <--SVI
  lines(slicefitdata[,1], slicefitdata[,7], col = slicecol[i]);
}
plot.new() # new plot for legend only
legend("center", rev(c("0.03","0.07","0.11", "0.18", "0.36", "0.61", "1.1", "1.6")), col=rev(slicecol), pch="o", 
       title = "T in years");

summederror
summedcrossedness

#PLOTTING DENSITY FOR BUTTERFLY
for (i in 1:numslices){
  slicesize <- 100;
  slicefitdata <- array(NA, dim=c(slicesize,7));							
  slicefitdata[,1] <- seq(from = -5, to = 3, length.out = slicesize);			#log moneyness k
  slicefitdata[,2] <- rep(surfacesliceparams[i,1], times = slicesize);		#ATM market total variance (proxy for time)
  slicefitdata[,3] <- rep(surfacesliceparams[i,2], times = slicesize);
  slicefitdata[,4] <- rep(surfacesliceparams[i,3], times = slicesize);
  slicefitdata[,5] <- rep(surfacesliceparams[i,4], times = slicesize);
  slicefitdata[,6] <- rep(surfacesliceparams[i,5], times = slicesize);
  slicefitdata[,7] <- apply(slicefitdata[,1:6], 1, bfdensity);
  if(i==1){
    plot(slicefitdata[,1], slicefitdata[,7], col = slicecol[i], type="l",xlim=c(-1.0, 1.0), ylim=c(-0.1,2.0), main="SSVI fit to volatility surface",
         xlab="Log moneyness k", ylab="Density");
    abline(a=0, b=0);
  }else{
    lines(slicefitdata[,1], slicefitdata[,7], col = slicecol[i]);
  }
}

plot.new() # new plot for legend only
legend("center", rev(c("0.03","0.07","0.11", "0.18", "0.36", "0.61", "1.1", "1.6")), col=rev(slicecol), pch="o", title = "T");



#5. OPTIMAL CALIBRATION: SVI FIT, PENALIZING FOR CROSSINGS AND BUTTERFLY ARBITRAGE	
#5a. FITTING SLICES WITH RAW SVI PARAMETRIZATION USING SIMULATED ANNEALING
optsliceparamsSA <- array(NA, dim=c(5,numslices)); #optimal slice parameters using Simulated Annealing method
summederror <- 0; # SVI price error
summedcrossedness <- 0; #crossed-ness
for (i in 1:numslices){
  
  t = DF_Flist$ttm[i]
  dataset_slice_SA = dataset.unique[which(dataset.unique$ttm==t), ];
  slicefit_SA <- NULL;
  startlist <- c(a = surfacesliceparams[i,1], b = surfacesliceparams[i,2], rho = surfacesliceparams[i,3], m = surfacesliceparams[i,4], sigma = surfacesliceparams[i,5]);
  if(i==1){
    slicefit_SA <- GenSA(par=startlist, fn=wrawpriceSA_bf , lower = c(a=-100,b=0,rho=-1,m=-100,sigma=0), upper = c(a=100,b=100,rho=1,m=100,sigma=100),
                         control=list(maxit=500, temperature=5e6, visiting.param=2.7, acceptance.param=-5), datavec=dataset_slice_SA);
  }else{
    slicefit_SA <- GenSA(par=startlist, fn=wrawpriceSA_x_bf , lower = c(a=-100,b=0,rho=-1,m=-100,sigma=0), upper = c(a=100,b=100,rho=1,m=100,sigma=100),
                         control=list(maxit=500, temperature=5e6, visiting.param=2.7, acceptance.param=-5), datavec=dataset_slice_SA, lastsliceparams=optsliceparamsSA[,i-1]);
  }
  if(is.null(slicefit_SA)){
    optsliceparamsSA[,i] <- rep(NA, times = 5);
  } else {
    optsliceparamsSA[,i] <- slicefit_SA$par;
  }
  summederror <- summederror + wrawpriceSA(optsliceparamsSA[,i], dataset_slice_SA);
  slice_k = log(dataset_slice_SA$strike/dataset_slice_SA$forward);	
  slice_w = ((dataset_slice_SA$ivol)^2)*t
  if(i==1){
    plot(slice_k, slice_w, col = slicecol[i], xlim=c(-0.2,0.15), 
         ylim = c(0,0.08), type="p", xlab = "Log moneyness", ylab = "Total implied variance" )
  } else{
    points(slice_k, slice_w, col = slicecol[i], type="p", xlab = "Log moneyness", ylab = "Total implied variance" )
    crossedness_arg = as.data.frame(rbind(optsliceparamsSA[c(1,2,5,3,4),i-1], optsliceparamsSA[c(1,2,5,3,4),i]));
    colnames(crossedness_arg)<- c("a","b","sig","rho","m");
    summedcrossedness <- summedcrossedness + ComputeSviRoots(crossedness_arg)$crossedness;
  }
  slicesize <- 100;
  slicefitdata <- array(NA, dim=c(slicesize,7));							
  slicefitdata[,1] <- seq(from = -2, to = 1, length.out = slicesize);			#log moneyness
  slicefitdata[,2] <- rep(optsliceparamsSA[1,i], times = slicesize);			#a
  slicefitdata[,3] <- rep(optsliceparamsSA[2,i], times = slicesize);			#b
  slicefitdata[,4] <- rep(optsliceparamsSA[3,i], times = slicesize);      #sig
  slicefitdata[,5] <- rep(optsliceparamsSA[4,i], times = slicesize);      #rho
  slicefitdata[,6] <- rep(optsliceparamsSA[5,i], times = slicesize);      #m
  slicefitdata[,7] <- apply(slicefitdata[,1:6], 1, wraw2);                 #w(k,a,b,...)
  lines(slicefitdata[,1], slicefitdata[,7], col = slicecol[i]);
}
summederror;
summedcrossedness;
optsliceparamsSA->optsliceparamsSAraw_x_bf; #optslicepramSA- raw - cross n butterfly penalty

#5b. PLOT 2D OF THE RAW FITS
for (i in 1:numslices){
  
  t = DF_Flist$ttm[i]
  dataset_slice_SA = dataset.unique[which(dataset.unique$ttm==t), ];
  slicefit_SA <- NULL;
  slice_k = log(dataset_slice_SA$strike/dataset_slice_SA$forward);	
  slice_w = ((dataset_slice_SA$ivol)^2)*t
  if(i==1){
    plot(slice_k, slice_w, col = slicecol[i], xlim=c(-0.5,0.4), ylim = c(0,0.08), 
         type="p", xlab = "Log moneyness k", ylab = "Total implied variance w(k,T)", 
         main="Raw slice fit (no static arbitrage)" )
  } else{
    points(slice_k, slice_w, col = slicecol[i], type="p");
  }
  slicesize <- 100;
  slicefitdata <- array(NA, dim=c(slicesize,7));							
  slicefitdata[,1] <- seq(from = -2, to = 1, length.out = slicesize);			#log moneyness
  slicefitdata[,2] <- rep(optsliceparamsSAraw_x_bf[1,i], times = slicesize);	#ATM market total variance (proxy for time)
  slicefitdata[,3] <- rep(optsliceparamsSAraw_x_bf[2,i], times = slicesize);
  slicefitdata[,4] <- rep(optsliceparamsSAraw_x_bf[3,i], times = slicesize);
  slicefitdata[,5] <- rep(optsliceparamsSAraw_x_bf[4,i], times = slicesize);
  slicefitdata[,6] <- rep(optsliceparamsSAraw_x_bf[5,i], times = slicesize);
  slicefitdata[,7] <- apply(slicefitdata[,1:6], 1, wraw2);
  lines(slicefitdata[,1], slicefitdata[,7], col = slicecol[i]);
}
plot.new() # new plot for legend only
legend("center", rev(c("0.03","0.07","0.11", "0.18", "0.36", "0.61", "1.1", "1.6")), col=rev(slicecol), pch="o", title = "T");

#5c. PLOTTING DENSITY FOR BUTTERFLY
for (i in 1:numslices){
  slicesize <- 500;
  slicefitdata <- array(NA, dim=c(slicesize,7));							
  slicefitdata[,1] <- seq(from = -5, to = 3, length.out = slicesize);			#log moneyness k
  slicefitdata[,2] <- rep(optsliceparamsSAraw_x_bf[1,i], times = slicesize);			#a
  slicefitdata[,3] <- rep(optsliceparamsSAraw_x_bf[2,i], times = slicesize);			#b
  slicefitdata[,4] <- rep(optsliceparamsSAraw_x_bf[3,i], times = slicesize);
  slicefitdata[,5] <- rep(optsliceparamsSAraw_x_bf[4,i], times = slicesize);
  slicefitdata[,6] <- rep(optsliceparamsSAraw_x_bf[5,i], times = slicesize);
  slicefitdata[,7] <- apply(slicefitdata[,1:6], 1, bfdensity);
  if(i==1){
    plot(slicefitdata[,1], slicefitdata[,7], col = slicecol[i], type="l", xlim=c(-3, 2.5), ylim=c(-0.1,3.5), 
         main="Raw slice fit (no static arbitrage)", xlab="Log moneyness k", ylab="Density");
    abline(a=0, b=0);
  }else{
    lines(slicefitdata[,1], slicefitdata[,7], col = slicecol[i]);
  }
}
plot.new() # new plot for legend only
legend("center", rev(c("0.03","0.07","0.11", "0.18", "0.36", "0.61", "1.1", "1.6")), col=rev(slicecol), 
       pch="o", title = "T")

optsliceparamsSAraw_x_bf<-t(optsliceparamsSAraw_x_bf) #transpose parameters matrix
colnames(optsliceparamsSAraw_x_bf)<-c("a","b",'rho','m', "sig")


