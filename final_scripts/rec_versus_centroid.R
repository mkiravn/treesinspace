

recursive_fraction <- function(k,j){
  result <- (k-1)/(k+j-2)
  return(result)
}

centroid_fraction <- function(k,j){
  result <- (k-1)*(k)/((k+j+1)*(k+j))
  return(result)
}


recursive_fraction(6,1)
centroid_fraction(6,1)

ratios <-  c(1:1000)/100


plot(x=ratios,y=recursive_fraction(ratios*100,100),type="l" ,xlab="k/j", ylab="fraction of comparisons within k")
points(x=ratios,y=centroid_fraction(ratios*100,100),col="red",type="l")
legend('topleft', legend=c("Recursive", "Centroid"),
       col=c("black", "red"), lty=1)

