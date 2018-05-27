################################################################################
#
# R code for the analysis of growth rate of rats by hierachical normal model.
# Follow the instructions given in the code to run requiored analyses.
#
################################################################################

# ---- Import data.csv file into R ----
data <- read.csv("data.csv")
# ---- Import data.csv file into R ----

# ---- After importing your data run the following codes from line 13 to line 87 ----
Y = as.matrix(data) # The y values are the component named y.

x = c(8.0, 15.0, 22.0, 29.0, 36.0)
xbar = mean(x)
N= nrow(Y)
t = ncol(Y)
dataList = list(x = x, xbar = xbar, N = N, t = t, Y= Y)

modelString = "
model
{
  for( i in 1 : N ) {
    for( j in 1 : t ) {
      Y[i , j] ~ dnorm(mu[i , j],tau.c)
      mu[i , j] <- alpha[i] + beta[i] * (x[j] - xbar)
    }
    alpha[i] ~ dnorm(alpha.c,alpha.tau)
    beta[i] ~ dnorm(beta.c,beta.tau)
  }
  tau.c ~ dgamma(0.001,0.001)
  sigma <- 1 / sqrt(tau.c)
  alpha.c ~ dnorm(0.0,1.0E-6)	   
  alpha.tau ~ dgamma(0.001,0.001)
  beta.c ~ dnorm(0.0,1.0E-6)
  beta.tau ~ dgamma(0.001,0.001)
  alpha0 <- alpha.c + xbar * beta.c	
}
"

writeLines( modelString , con="TEMPmodel.txt" ) # write to file

thetaInit = list(alpha = c(250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 
               250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250, 250),
     beta  = c(6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 
               6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6),			
     alpha.c = 150, beta.c = 10, 
     tau.c = 1, alpha.tau = 1, beta.tau = 1)

library(rjags)

jagsModel = jags.model( file="TEMPmodel.txt" , 
                        # the name of the file in which the model 
                        # specification is stored
                        data=dataList ,        
                        # the list of data
                        inits=thetaInit ,      
                        # the list of initial values
                        # to let JAGS to create its own 
                        # initial values for the chains, 
                        # simply omit this argument entirely
                        n.chains=3 ,           
                        # the number of chains to be generated
                        n.adapt=500            
                        # the number of steps to take for adapting 
                        # (or tuning) the samplers
)

update( jagsModel ,   # tell the name of the object 
        # that include the model to JAGS
        n.iter=500    # specify the length of the burn-in period
)

codaSamples = coda.samples( jagsModel ,                 
                            # previously created JAGS model object
                            variable.names=c("alpha0","beta.c","sigma") , 
                            # specify which parameters will have 
                            # their values recorded during the 
                            # MCMC walk
                            n.iter=3334                 
                            # specify the number of iterations for 
                            # each chain
)

source("DBDA2E-utilities.R")
# ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

# --- Run the following line to get diagnostic plots for overall mean growth rate ---
diagMCMC( codaObject=codaSamples , parName="alpha0" )

# --- Run the following line to get diagnostic plots for mean of unit growth ---
diagMCMC( codaObject=codaSamples , parName="beta.c" )

# --- Run the following line to get diagnostic plots for overall precision ---
diagMCMC( codaObject=codaSamples , parName="sigma" )


# --- Run the following plotPost() function to get posterior histogram and posterior estimates for overall mean growth rate ---
plotPost( codaSamples[,"alpha0"] , # the element of the posterior 
          # samples to be plotted
          main="theta" ,          # main title    
          xlab=bquote(alpha0)      # x-axis label
)

# --- Run the following plotPost() function to get posterior histogram and posterior estimates for mean of unit growth ---
plotPost( codaSamples[,"beta.c"] , # the element of the posterior 
          # samples to be plotted
          main="theta" ,          # main title    
          xlab=bquote(beta.c)      # x-axis label
)

# --- Run the following plotPost() function to get posterior histogram and posterior estimates for overall precision ---
plotPost( codaSamples[,"sigma"] , # the element of the posterior 
          # samples to be plotted
          main="theta" ,          # main title    
          xlab=bquote(sigma)      # x-axis label
)


