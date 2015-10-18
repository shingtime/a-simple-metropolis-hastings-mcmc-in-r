## generate the true model's slope
trueA <- 5
## generate the true model's intercept
trueB <- 0
## generate the true model's error term's standard deviation
trueSd <- 10

## generate the number of observations
sampleSize <- 31

# create independent x-values 
## generate x from -15 to 15
x <- (-(sampleSize-1)/2):((sampleSize-1)/2)
# create dependent values according to ax + b + N(0,sd)
## generate the true response variable for the linear relation
y <-  trueA * x + trueB + rnorm(n=sampleSize,mean=0,sd=trueSd)
## draw the scatter plot of relation between x and y
plot(x,y, main="Test Data")

source('all_functions_source.r')
# Example: plot the likelihood profile of the slope a
## when changing slope, the sum of log likelihood value changes
slopevalues <- function(x){return(likelihood(c(x, trueB, trueSd)))}
## vectorized caculating the log likelihood value for y when changing slope from 3 to 7 by 0.5.
slopelikelihoods <- lapply(seq(3, 7, by=.05), slopevalues )

## observe the relation between slope and the sum of log likelihood value
## find that when slope = 5, the log likelihood has the largest value, which is true since the true
## y is generating with slope 5.

plot (seq(3, 7, by=.05), slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")

# Prior distribution
## for  slope,intercept and standard deviation
# return the posterior loglikehood density value for slope,intercept and sd when observing y.


######## Metropolis algorithm ################
## choose stard value of slope, intercept and sd.
startvalue = c(4,0,10)
## generate all the possible sets of slope, intercept and sd using algorithms.
chain = run_metropolis_MCMC(startvalue, 10000)
## drop the first 5000 values since the first steps may be biased.
burnIn = 5000

## return the acceptance rate: how many cases do we accept they are from our posterior
## duplicated returns bool values
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

### Summary: #######################
## adjust plot areas: plot 2*3 plots in one window

par(mfrow = c(2,3))
## call the function with different arguments to see the results
plot_posterior(chain,burnIn,"slope",trueA,trueB,trueSd,TRUE)
plot_posterior(chain,burnIn,"intercept",trueA,trueB,trueSd,TRUE)
plot_posterior(chain,burnIn,"sd",trueA,trueB,trueSd,TRUE)
plot_posterior(chain,burnIn,"slope",trueA,trueB,trueSd,FALSE)
plot_posterior(chain,burnIn,"intercept",trueA,trueB,trueSd,FALSE)
plot_posterior(chain,burnIn,"sd",trueA,trueB,trueSd,FALSE)

# for comparison:
## using our y values to fit the linear regression in order to compare our true parameters.
summary(lm(y~x))


### compare_outcomes

compare_outcomes = function(iterations){
  ## generate a list to store 10 times results
  result = list()
  ## specify how many obs will be discard
  burnIn = 5000 
  ## loop 10 times
  for (i in 1:10){
    ## for each iteration, generate startvalues which follows the same prior distribution as specified
    ## before.
    startvalue = c(runif(1,0,10),rnorm(1,0,5),runif(1,0,30))
    ## generate the chain
    chain = run_metropolis_MCMC(startvalue, 10000)
    ## get all the values for a
    a = chain[-(1:burnIn),1]
    #a = chain[,1]
    ## store mean and sd into a list
    result[[i]] = c(mean = mean(a),sd = sd(a)) 
  }
  ## return 10 times result
  result
}

## test cases

compare_outcomes(1000)
compare_outcomes(10000)
compare_outcomes(100000)


