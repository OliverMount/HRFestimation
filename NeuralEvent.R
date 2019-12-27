NeuralEvent<-function(N,P,ty)
{
  # March 15, 2019
  # This functions returns the Neural vector depends on the type  

Imp<-array(0,dim=c(N))   # Initialize neural events vector

if (size(P)[2]==2)     # If there are negative events
{
  Peve<-sort(P[[1]])    # Positive events
  Neve<-sort(P[[2]])    # Record Negative events
} else {
    Peve<-sort(P)
  }   

switch(ty,
       "P" =
       {
         Imp[Peve]<- 1    
       },
       
       "PN" =
       {
         Imp[Peve]<- 1      # Positive impulse
         Imp[Neve]= -1      # Negative Impulse
       },
       
       "PM" =
       {
         Imp[Peve]<- 1      # Positive impulse
       },
       
       "PNM" =
       {
         Imp[Peve]<- 1      # Positive impulse
         Imp[Neve]= -1      # Negative Impulse
       }
       
)
return(Imp)
}