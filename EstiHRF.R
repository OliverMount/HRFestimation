EstiHRF<-function(sig,Imp,IBa,TR)
{
# This function returns the estimated HRF 

LBB<-dim(IBa)[1]  # length of the Basis functions

## Finding optimal delay and HRF
kmax=15     # maximum delay of the hemodynamic response  to search for
err<-array(NA,dim=c(kmax))

for (k in 1:kmax)
{  
  # sig= C*h +n  This is the convolution model assumed
  Imp1<-shift(Imp,n=(k-1),fill=0,type="lead")    # shift the neural impulses ( that occur before the Hemodynamics)
  C<-convolmtx(Imp1,LBB)        # Convolution matrix of impulse response
  X <- C%*%IBa         # Design matix (without intercept)
  Xi <- cbind(X,dim=c(ones(dim(X)[1],1))) # Design matix (with intercept)
  
  D<-Xi[1:(N),]   
  h<-pinv(D) %*% sig
   #ModeledSig= D%*%h  # Estimated signal  (with intercept)
  ModeledSig= X[1:(N),]%*%h[1:dim(X)[2]]  # Estimated signal (without intercept)
  err[k]= sum((sig-ModeledSig)^2)        # squared error (of noise!)
  # plot(sig,type="l",col="blue")
  # lines(ModeledSig,col="red")
}

# Kopt
Kopt<-which(err==min(err)) # in samples

# Modeled signal
Imp1<-shift(Imp,n=(Kopt-1),fill=0,type="lead")    # shift the neural impulses ( that occur before the Hemodynamics)
C<-convolmtx(Imp1,LBB)        # Convolution matrix of impulse response
X <- C%*%IBa         # Design matix (without intercept)
Xi <- cbind(X,dim=c(ones(dim(X)[1],1))) # Design matix (with intercept)
D<-Xi[1:(N),]   
h<-pinv(D) %*% sig
ModeledSig= X[1:(N),]%*%h[1:dim(X)[2]]

# Estimated HRF
HRF<-(IBa%*%h[1:dim(X)[2]])  

 # plot(sig,type="l",col="blue",ylim=c(-3.5,3.5))
 # lines(ModeledSig,col="red")

Kopt<- (Kopt*TR) # in sec  
return(list(HRF=HRF,ModeledSig=ModeledSig,Kopt=Kopt))

}