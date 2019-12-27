convolmtx <- function(x,M)
{ 
# CONVOLMTX(X,M) Returns the Length(x)+M-1 X M convolution matrix 
# x- Input data sequence -  in row wise input.
# M- length of the filter
  
r<-array(0,dim=c(length(x)+M-1,M))

for (i in 1:M)
{
r[,i]=c(rep(0,i-1),x,rep(0,M-i))  #  column wise implementation 
} 
return(r)
}