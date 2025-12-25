

#1. BS: price and iv
#2. Function for fiting with ssvi price: initial guess
#3. Function for fitting SVI price: optimal calibration (wrawprice)


#1a. BS: price. price of call option 

BScallprice<- function(ttm, strike, sigma, DF, forward) #BS price.callprice at t=0, for a option at T
{
  k <- log(strike/forward); 
  sig <- sigma*sqrt(ttm); 
  d1 <- -k/sig+sig/2; #it is the same as d1 in BS formula
  d2 <- d1 - sig;
  callprice = DF*(forward*pnorm(d1) - strike*pnorm(d2));										#Black-Scholes price of a callprice
  callprice
}
#try: BScallprice( ttm=0.05, strike=2350, sigma=0.40040, DF=0.9985929, forward=2850.018)
#BScallprice( ttm=1/365, strike=80, sigma=0.576, DF=exp(-0.0025/365), forward=83.11*exp(0.0025/365)) price =3.23

#1b. IV COMPUTE
# This function now works with vectors of strikes and option values. GET IV
BSImpliedVolCall <- function(ttm, strike, DF, forward, C) #C is then: BSprice = C --> IV
{
  nK <- length(strike);
  sigmaL <- rep(1e-10,nK); #set value for iv (sigma). first, set min value = 1e-10
  CL <- BScallprice(ttm, strike, sigmaL, DF, forward); #compute the corresponding BSprice for that iv
  sigmaH <- rep(10,nK);
  CH <- BScallprice(ttm, strike, sigmaH, DF, forward);
  while (mean(sigmaH - sigmaL) > 1e-10)
  {
    sigma <- (sigmaL + sigmaH)/2;
    CM <- BScallprice(ttm, strike, sigma, DF, forward);#Bs price, mid
    CL <- CL + (CM < C)*(CM-CL);
    sigmaL <- sigmaL + (CM < C)*(sigma-sigmaL);
    CH <- CH + (CM >= C)*(CM-CH);
    sigmaH <- sigmaH + (CM >= C)*(sigma-sigmaH);
  }
  return(sigma);
}

#2. fitting SSVI to prices: functions: putprice of surface
SSVIprice <- function(ttm, strike, DF, forward, w0, eta, rho ){ 
 
  k = log(strike/forward); #x =^= log moneyness, w0 =^= ATM total implied var
  phi= eta/sqrt(w0*(w0+1)) #phi follows complex power law form 
  w = w0/2* (1+rho*phi*k + sqrt((phi*k + rho)^2 +(1-rho^2) ) );	#total var SSVI form. 
  sig <- sqrt(w); 																				
  d1 <- -k/sig+sig/2; #it is the same as d1 in BS formula
  d2 <- d1 - sig;
  callprice = DF*(forward*pnorm(d1) - strike*pnorm(d2));										#Black-Scholes price of a callprice
  callprice;															
}

SSVIprice2 <- function(invec){
  ttm =invec[1]; strike = invec[2]; DF = invec[3]; forward = invec[4]; w0 = invec[5]; eta = invec[6]; rho = invec[7];	
  SSVIprice(ttm, strike, DF, forward, w0, eta, rho);
}
SSVIpriceSA <- function(paramvec, datavec){ 
  eta = paramvec[1]; rho = paramvec[2]; 
  datavec_length = dim(datavec)[1];
  totalerror = 0; 
  for (i in 1:datavec_length){
  tem1= c(as.numeric(datavec[i,1]),as.numeric(datavec[i,2]),as.numeric(datavec[i,3]),as.numeric(datavec[i,4]),
           as.numeric(datavec[i,5]), eta, rho)
    
  totalerror = totalerror + (SSVIprice2(tem1)-datavec[i,6])^2/datavec[i,6];
  }
  totalerror = totalerror + max(0,eta*(1+abs(rho))-2); #one penalty because we just need ssvi satisify one condition
  totalerror;
}

#SSVI function (total implied volatility surface based on 2 params)
SSVI <- function(k, w0, eta, rho){
  phi= eta/sqrt(w0*(w0+1)) #phi follows complex power law form 
  w0/2* (1+rho*phi*k + sqrt((phi*k + rho)^2 +(1-rho^2) ) );	#total var SSVI form.
}
SSVI2 <- function(invec){
  k = invec[1]; w0=invec[2];eta=invec[3];rho=invec[4];
  SSVI(k,w0,eta,rho)
}
#2b. Parameter convert: SSVI--> SVI-JW --> Raw
#first define parameters from SSVI to RAW and raw parametrization itself
SSVItoJW <- function(eta, rho, w0, ttm){
  v <- w0/ttm;
  psi <- 0.5*rho*sqrt(w0)*eta/sqrt(w0)/sqrt(1+w0);
  p <- 0.5*sqrt(w0)*eta/sqrt(w0)/sqrt(1+w0)*(1-rho);
  cc <- 0.5*sqrt(w0)*eta/sqrt(w0)/sqrt(1+w0)*(1+rho);
  vhat <- w0/ttm*(1-rho^2);
  output <- c(v,psi,p,cc,vhat);
  names(output) <- c("v", "psi", "p", "c", "vhat");
  output;
}
JWtoRAW <- function(invec){
  v <- invec[1]; psi <- invec[2]; p <- invec[3]; cc <- invec[4]; vhat <- invec[5]; ttm <- invec[6];
  w <- v*ttm;
  b <- 0.5*sqrt(w)*(cc+p);
  rho <- 1-p*sqrt(w)/b;
  bbeta <- rho-2*psi*sqrt(w)/b;
  aalpha <- sign(bbeta)*sqrt(1/bbeta^2-1);
  m <- (v-vhat)*ttm/b/(-rho+sign(aalpha)*sqrt(1+aalpha^2)-aalpha*sqrt(1-rho^2));
  sigma <- aalpha*m;
  a <- vhat*ttm-b*sigma*sqrt(1-rho^2);					
  if(m == 0){sigma <- (v*ttm-a)/b};							#well but this uses so far undefined a which itself requires sigma
  output <- c(a,b,rho,m,sigma);
  names(output) <- c("a", "b", "rho", "m", "sigma");
  output;
}

#3. wrawprice: 4 versions: no penalty for arb, add penalty for calendar arb, bf arb, both arb
#3a: wrawpriceSA: total error. 
#3b: wrawpriceSA_x: total error + penalty.cross
#3c: wrawpriceSA_bf: total error + penalty.bf
#3d: wrawpriceSA_x_bf: total error + penalty.cross+ penalty.bf

#3a:wrawpriceSA: total error.
wraw <- function(k,a,b,rho,m,sigma){
  #a = chi[1];b = chi[2]; rho = chi[3]; m = chi[4]; sigma = chi[5];
  wout = a+b*(rho*(k-m)+sqrt((k-m)^2+sigma^2));
  wout;
}
wraw2 <- function(chi){					#chi is the parameter vector
  k=chi[1]; a = chi[2];b = chi[3]; rho = chi[4]; m = chi[5]; sigma = chi[6];
  wout = wraw(k,a,b,rho,m,sigma);
  wout;
}
wrawprice <- function( ttm, strike, DF, forward, a, b, rho, m, sigma){
  k <- log(strike/forward);
  #w = a+b*(rho*(k-m)+sqrt((k-m)^2+sigma^2));
  w = wraw(k, a, b, rho, m, sigma);
  if(w < 0){w=0};
  sig = sqrt(w) 
  d1 <- -k/sig+sig/2; #it is the same as d1 in BS formula
  d2 <- d1 - sig;
  callprice = DF*(forward*pnorm(d1) - strike*pnorm(d2));										#Black-Scholes price of a call
  callprice
}
wrawprice2 <- function(invec){
  ttm = invec[1]; strike = invec[2]; DF = invec[3]; forward = invec[4]; a = invec[5];b = invec[6]; rho = invec[7]; m = invec[8]; sigma = invec[9];
  callprice = wrawprice(ttm, strike, DF, forward, a, b, rho, m, sigma);						#Black-Scholes price of a put
  callprice;
}
wrawpriceSA <- function(paramvec, datavec){
  a = paramvec[1]; b = paramvec[2]; rho = paramvec[3]; m = paramvec[4]; sigma = paramvec[5];
  datavec_length = dim(datavec)[1];
  totalerror = 0; 

  for (i in 1:datavec_length){
    tem2= c(as.numeric(datavec[i,1]),as.numeric(datavec[i,2]),as.numeric(datavec[i,3]),as.numeric(datavec[i,4]),
            a, b, rho, m, sigma)
    totalerror = totalerror + (wrawprice2(tem2)-datavec[i,6])^2/datavec[i,6];
  }
  #totalerror = totalerror #+penalty
  totalerror;
}


#3b: wrawpriceSA_x: total error + penalty.cross
wrawpriceSA_x <- function(paramvec, datavec, lastsliceparams){
  a = paramvec[1]; b = paramvec[2]; rho = paramvec[3]; m = paramvec[4]; sigma = paramvec[5];
  datavec_length = dim(datavec)[1];
  totalerror = 0; 
  for (i in 1:datavec_length){
    tem3= c(as.numeric(datavec[i,1]),as.numeric(datavec[i,2]),as.numeric(datavec[i,3]),as.numeric(datavec[i,4]),
            a, b, rho, m, sigma)
    totalerror = totalerror + (wrawprice2(tem3)-datavec[i,6])^2/datavec[i,6];
  }
  crossedness_arg = as.data.frame(rbind(lastsliceparams[c(1,2,5,3,4)], c(a,b,sigma, rho, m)));
  colnames(crossedness_arg)<- c("a","b","sig","rho","m");
  penalty = ComputeSviRoots(crossedness_arg)$crossedness; 
  totalerror = totalerror + penalty#10*penalty;
  totalerror;
}

#3c: wrawpriceSA_bf: total error + penalty.bf

wrawp <- function(chi){						#first deriv
  k=chi[1]; a = chi[2];b = chi[3]; rho = chi[4]; m = chi[5]; sigma = chi[6];
  b*rho+b*(k-m)/sqrt((k-m)^2+sigma^2);
}
wrawpp <- function(chi){					#second deriv of wraw wrt k
  k=chi[1]; a = chi[2];b = chi[3]; rho = chi[4]; m = chi[5]; sigma = chi[6];
  b*sigma^2/((k-m)^2+sigma^2)^(3/2);
}
#density: if negative butterfly arbitrage is possible
bfdensity <- function(chi){
  k=chi[1]; a = chi[2];b = chi[3]; rho = chi[4]; m = chi[5]; sigma = chi[6];
  w = wraw2(chi);
  wp = wrawp(chi);
  wpp = wrawpp(chi);
  wout = (1-k*wp/2/w)^2-wp^2/4*(1/w+1/4)+wpp/2;
  wout;
}
wrawpriceSA_bf <- function(paramvec, datavec){
  a = paramvec[1]; b = paramvec[2]; rho = paramvec[3]; m = paramvec[4]; sigma = paramvec[5];
  datavec_length = dim(datavec)[1];
  totalerror = 0; 
  for (i in 1:datavec_length){
    tem4= c(as.numeric(datavec[i,1]),as.numeric(datavec[i,2]),as.numeric(datavec[i,3]),as.numeric(datavec[i,4]),
            a, b, rho, m, sigma)
    totalerror = totalerror + (wrawprice2(tem4)-datavec[i,6])^2/datavec[i,6];
  }
  bfpenalty = 0;
  for (i in seq(from = -5, to = 3.0, length.out=50)){
    bfpenalty = bfpenalty + abs(min(0, bfdensity(c(i, paramvec))));
  }
  totalerror = totalerror + bfpenalty #0.1*bfpenalty;
  totalerror;
}

#3d: wrawpriceSA_x_bf: total error + penalty.cross+ penalty.bf
wrawpriceSA_x_bf <- function(paramvec, datavec, lastsliceparams){
  a = paramvec[1]; b = paramvec[2]; rho = paramvec[3]; m = paramvec[4]; sigma = paramvec[5];
  datavec_length = dim(datavec)[1];
  totalerror = 0; 
  for (i in 1:datavec_length){
    tem5= c(as.numeric(datavec[i,1]),as.numeric(datavec[i,2]),as.numeric(datavec[i,3]),as.numeric(datavec[i,4]),
            a, b, rho, m, sigma)
    totalerror = totalerror + (wrawprice2(tem5)-datavec[i,6])^2/datavec[i,6];
  }
  crossedness_arg = as.data.frame(rbind(lastsliceparams[c(1,2,5,3,4)], c(a,b,sigma, rho, m)));
  colnames(crossedness_arg)<- c("a","b","sig","rho","m");
  penalty = ComputeSviRoots(crossedness_arg)$crossedness; 
  bfpenalty = 0;
  for (i in seq(from = -5, to = 3.0, length.out=50)){
    bfpenalty = bfpenalty + abs(min(0, bfdensity(c(i, paramvec))));
  }
  totalerror = totalerror + 1*penalty + 1*bfpenalty;#10*penalty
  totalerror;
}
