Give_HRF_para<-function(sig,TR,MII,ty,LB,Thres=1)
{
  
  # This functions returns the point process (pseudo events) and 
  # Shape of HRF along with its parameters (FWHM, TTP, Delay,RH,TA,T2T)
  # The delalyed is chosen such the error noise variance is the least
  
  # Input parameters (list in the following order)
  # sig < - Voxel time series
  # LB <-  length of the BOLD HRF (length of the HRF)
  # TR <-   Repetition time in fMRI (sampling time)
  # MII   - Multiple Index Identifier
  #         A value (in TRs) below which multiple peaks are less likely (defaults to 12 TRs)        
  #         (Please see function Give_Events.R for more details about this parameter)  
  
  # output parameters (list in the following order)
  
  # 1. para[1]  -- Pseudo events index 
  # 2. para[2]  -- HRF
  # 3. para[3]  -- optimal delay
  # 4. para[4]  -- Full width half maximum (FWHM)
  # 5. para[5]  -- Time to peak            (TTP)
  # 6. para[6] --  Response height         (RH)
  # 7. para[7]  -- First zero crossing in HRF for verifying the TTP results
  # 8. para[8]  -- Scaling factor
  # 9. para[9]  -- Modeled signal
  
  library(BHMSMAfMRI) # For handling fMRI dataset
  library(signal)  # for Windows and filters
  library(data.table)  # For shifting
  library(DescTools)  # Fisher Z score
  library(mosaic)   # Z- score
  library(SuppDists)  # For Gamma distribution test for Rho
  library(Matrix)  # For quadratic form
  library(devtools)
  library(pracma)  # For pinv, fprintf, size
  library(fields)  # For yline, xline
  library(FIACH)  # For informed basis functions and spm stuffs
  
  ### Please change the path for other functions here
  source("Give_Events.R")
  source("IBasis.R")
  source("NeuralEvent.R")
  source("EstiHRF.R")  
  source("find_peaks.R")
  
  ####################   Finding Pseudo events (positive+negative+multiple peaks by default)  #################### 
  P<- Give_Events(sig,ty,MII,Thres)   # see the function Give_Events for more details about the type
  sig<-zscore(sig)
  ####################  Estimation of HRF by spike triggered average method #################### 
  
  N<-length(sig)
  
  ##### HRF estimation  starts here #### 
  
  # Get an orthogonal informed basis (given TR) 
  IBa<-IBasis(TR,orth=T,LB)
  LBB<-dim(IBa)[1]  # length of the Basis functions
  Imp<-NeuralEvent(N,P,ty)  # Neural event vector (P contains the timings)
  HRF_EstiSig <- EstiHRF(sig,Imp,IBa,TR)  # HRF and Estimated signal returned
  HRF<- HRF_EstiSig$HRF
  ModelSig<- HRF_EstiSig$ModeledSig
  Kopt<- HRF_EstiSig$Kopt
  
  if (!length(Kopt))
  {Kopt<-0}  
  
  ####################  Estimation of HRF parameters #################### 
  
  ## FWHM
  h_ma<-0.5*max(HRF)
  temp<-which(HRF>=h_ma)
  FWHM<-(tail(temp,n=1)-temp[1])  # in samples
  FWHM<- (FWHM * TR) #  in sec
  
  # Response Height
  RH<-max(HRF)
  
  # Time to peak
  ttp<-TTP<-which(HRF==RH)  # in samples
  TTP<- (TTP*TR)   #  in sec
  
  # Other parameters
  
  # Trough Amplitude (TA)
  TA<-min(HRF[ttp:LBB])  # Find the min between peak and end
  
  # Time to Trough (T2T)
  T2T<-which(HRF==TA)  # in samples
  T2T<-(T2T*TR)
  
  # Peak to Trough (P2T)
  P2T<-(T2T-TTP)    #  in sec
  
  
  ####################  return the desired parameters to main program   #################### 
  # Make a list of the desired parameters 
  # In order: Events, HRF,delay, FWHM, TTP, RH, Modeled signal, TA, T2T, P2T
   
  return(list(P,HRF=HRF,Kopt=Kopt,FWHM=FWHM,TTP=TTP,RH=RH,ModeledSig=ModelSig,TA=TA,T2T=T2T,P2T=P2T))
  
}