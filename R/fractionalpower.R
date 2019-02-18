fractionalpower <- function(t,qt){
  x <- log10(t)
  y <- log10(qt)
  
  fit <- lm(y~x)
  print(summary(fit))
  
  rhs <- function(x,b0,b1){
    log10(b0) + b1*x
  }
   fit1 <- nlsLM(y~rhs(x,a,b), start=list(a=1,b=1))
   print("parameters")
   print(summary(fit1))
}
