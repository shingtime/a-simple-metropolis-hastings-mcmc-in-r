
likelihood <- function(param){  # input is the parameter vector
  a = param[1]   # true slope
  b = param[2]   # true intercept
  sd = param[3]  # true standard deviation
  
  pred = a*x + b # the expectation of y 
  singlelikelihoods = dnorm(y, mean = pred, sd = sd, log = T)  # return the corresponding log density if y has mean pred and variance sd^2
  sumll = sum(singlelikelihoods) ## since we take the log of PDF value, the joint density value transfer from product to sum
  return(sumll)   # return the sum of log likelihood values for true y
}
# Prior distribution
## for  slope,intercept and standard deviation
prior <- function(param){
  # quantile to generate slope parameter
  a = param[1]
  # quantile to generate intercept parameter
  b = param[2]
  # quantile  to generate standard deviation parameter
  sd = param[3]
  # specify slope follows a uniform distribution in(0,10) with quantile a , and return the log PDF value for slope
  aprior = dunif(a, min=0, max=10, log = T) 
  # specify intercept follows a normal distribution with mean = 0 and sd = 5. return the log PDF value
  bprior = dnorm(b, sd = 5, log = T)
  # specify sd follows uniform distribution in (0,30) with quantile sd, return the log PDF value
  sdprior = dunif(sd, min=0, max=30, log = T)
  
  ## return the sum of each quantile's log likelihood value for all parameters' prior distribution
  return(aprior+bprior+sdprior)
}
# return the posterior loglikehood density value for slope,intercept and sd when observing y.
posterior <- function(param){
  # return the sum of log likelihood value 
  return (likelihood(param) + prior(param))
}
######## Metropolis algorithm ################
## function to generate new parameter value close to old parameter value
proposalfunction <- function(param){
  ## generate 3 potential parameter with different sd but given mean
  return(rnorm(3,mean = param, sd= c(0.1,0.5,0.3)))
}


## implement the metropolis algorithms, given startvalue of parameters and iteration times
run_metropolis_MCMC <- function(startvalue, iterations){
  ## geneate array to store the startvalue and the calculated values in the loop(may update or may keep the old value)
  chain = array(dim = c(iterations+1,3))
  ## store the start value at the first position
  chain[1,] = startvalue
  ## loop to implement metropolis algorithms
  for (i in 1:iterations){
    ## using the startvalue and jumping distribution to generate a proposed parameter value
    proposal = proposalfunction(chain[i,])
    ## calculated the acceptance ratio(since we take log, there should be minus instead of divide)
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    ## criterion if accept new value follows our given posterior for parameters
    if (runif(1) < probab){
      ## if the ratio is less than a random number in (0,1) with a uniform distribution, accept the proposed value
      chain[i+1,] = proposal
    }else{
      ## else, keep the old value
      chain[i+1,] = chain[i,]
    }
  }
  ## for given iteration times, return the array of storing all the possible values of slope,intercept and sd.
  return(chain)
}

## function for the plots in the summary part
##chain,burnIn,trueA,trueB and trueSd are parameters that are given previously
## choice_of_plot asks user to input which parameter's plot he wants to show. There are three 
## choices: "slope","intercept" and "sd".
## choice_of_hist is a boolean variable. If it equals to TRUE, the function plots the usual
## histgram with estimator and true value as vertical lines. If this variable equals to FALSE,
## the function plots the chain value with the corresponding true value.

plot_posterior = function(chain,burnIn,choice_of_plot,trueA,trueB,trueSd,choice_of_hist){
  if(choice_of_hist==TRUE){
    if (choice_of_plot=="slope"){
      hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line" )
      abline(v = mean(chain[-(1:burnIn),1]))  
      abline(v = trueA, col="red" )
    }
    else if (choice_of_plot=="intercept"){
      hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
      abline(v = mean(chain[-(1:burnIn),2]))  
      abline(v = trueB, col="red" )
    }
    
    else{
      hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
      abline(v = mean(chain[-(1:burnIn),3]) )
      abline(v = trueSd, col="red" )
      
    }
  }
  
  if(choice_of_hist==FALSE){
    if (choice_of_plot=="slope"){
      plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
      abline(h = trueA, col="red" )
    }
    else if(choice_of_plot=="intercept"){
      plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
      abline(h = trueB, col="red" )
    }
    else{
      plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
      abline(h = trueSd, col="red" )
    }
  } 
}


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

