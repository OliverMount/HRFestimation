Give_Events<-function(sig,type,MII=12,Thres=1)
{
  
  # Inputs
  # Sig  - Time series
  # type  - Type of events needed
  #         P - positive only events
  #         PN - Positive + negative events
  #         PM - Positive + multipeak events
  #         PNM - Positive + negative + multipeak events (Default)
  # MII   - Multiple Index Identifier
  #         A value (in TRs) below which multiple peaks are less likely (defaults to 12 TRs)          
  # Thres  - Threshold to identity the events (default one standard deviation)
  
  # Check for zero-valued voxels before Z scoring the values
  
  if (!(sum(sig)==0)) # To exclude voxels with zero values (in the brain mask)
  {sig<-zscore(sig)} # Normalize the signal intensities by zscores  
  
  # To find only postive events (This has to be found for all the four cases)
  xp<- which( sig >= Thres*sd(sig))                     # one sd. threshold
  Ind1=c(1,(which(diff(xp)!=1)+1),(length(xp)+1))                 # Find positive indces (location)
  P<-xp[Ind1]

  switch(type,
         
         "PN" = {  # To find only negative events and append to the postive events
                   xn<- which( sig <= -Thres*sd(sig))           # one sd. threshold
                   Ind2=c(1,(which(diff(xn)!=1)+1),(length(xn)+1))  # Find negative indces (location)
                   Ne<-xn[Ind2]
                   P<-list(P,Ne)   # Append the negative part to the positive part
                }, 
         
         "PM" = {  # To find multiple peaks on the positive side
           
                   MIP<- which(diff(Ind1)>= MII)      # Find interval on positive side ( intervals in xp)
                   LMIP<-length(MIP) 
  
                       for (Mi in 1:LMIP)
                       {
                         Lindex<-xp[Ind1[MIP[Mi]]:(Ind1[MIP[Mi]+1]-1)]
                         Lpeaks<-find_peaks(-sig[Lindex],1)   # Finding minimum
             
                         if (!isempty(Lindex[Lpeaks]))
                         {
                           P<-c(P,Lindex[Lpeaks]) # Append index of multiple events  with the postive part
                         }                                           
                       }
           
                }, 
         
         "PNM" = {  # To find multiple peaks on both positive and negative side
  
                 # Multiple peaks on positive side
                     
                     MIP<- which(diff(Ind1)>= MII)      # Find no. of intervals greater than MII samples
                      
                     if (!isempty(MIP))    # If there is a multiple postive peaks then proceeed
                     { LMIP<-length(MIP)
                     for (Mi in 1:LMIP)
                     {
                       Lindex<-xp[Ind1[MIP[Mi]]:(Ind1[MIP[Mi]+1]-1)]
                       Lpeaks<-find_peaks(-sig[Lindex],1)   # Finding minimum
                       
                       if (!isempty(Lindex[Lpeaks]))
                       {
                         P<-c(P,Lindex[Lpeaks]) # Append index of multiple events  with the positive part
                       }                                           
                     }
                     }
                # Only negative side   
                     
                     xn<- which( sig <= -Thres*sd(sig))           # one sd. threshold
                     Ind2=c(1,(which(diff(xn)!=1)+1),(length(xn)+1))  # Find negative indces (location)
                     Ne<-xn[Ind2]
                     
                # Multiple peaks on negative side
           
                     MIN<- which(diff(Ind2)>= MII)       # Find no. of intervals greater than MII samples
                     
                     if (!isempty(MIN))    # If there is no multiple postive peaks then proceeed
                     { LMIN<-length(MIN)
                     for (Mi in 1:LMIN)
                     {
                       Lindex<-xn[Ind2[MIN[Mi]]:(Ind2[MIN[Mi]+1]-1)]
                       Lpeaks<-find_peaks(sig[Lindex],1)     # Finding maximum (because it is on the negative side)
                       
                       if (!isempty(Lindex[Lpeaks]))
                       {
                         Ne<-c(Ne,Lindex[Lpeaks])  # Append index of multiple events  with the positive  and negative part
                       }
                     }
                     }
                     
                     
                     P<-list(P,Ne) # Accumulate the positive and negative parts separately
         } 
        
)

  return(P)    # P is the index set of events
   
}