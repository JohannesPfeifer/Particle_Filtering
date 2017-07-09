# PARTICLE FILTERING AND SMOOTHING EXAMPLE CODE

These example codes illustrate the methods used in Benjamin Born/Johannes Pfeifer (2014): "Policy Risk and the Business Cycle", Journal of Monetary Economics, 68, pp. 68-85.  
Feel free to modify and adapt the codes to your needs, but please be fair and acknowledge the source. We ourselves have profited from the 
particle filter implementation of Andreasen, Martin M. (2011): "Non-Linear DSGE Models and The Optimized Central Difference Particle Filter", 
Journal of Economic Dynamics and Contol, 35(10), pp. 1671-1695


The following files do:
 1. Simulate AR1-Stochastic volatility process
 2. Run Metropolis-Hastings algorithm on AR1-stochastic volatility model using bootstrap (SIR) particle filter for evaluating the likelihood
 3. Run Doucet et. al. particle smoother

## How to run
The main file to run is: run_filter_and_smoother_AR1.m

## References: 

The particle filter follows Arulampalam/Maskell/Gordon/Clapp (2002): "A Tutorial on Particle Filters for Online Nonlinear/Non-Gaussian Bayesian Tracking", 
	IEEE Transactions on Signal Processing, 50(2)  

The smoother is implemented according to Godsill/Doucet/West(2004) "Monte Carlo smoothing for nonlinear time series", 
	Journal of the American Statistical Association, 2004, 99, 156-168

In case of questions or bugs, email us at jpfeifer@uni-koeln.de or born@uni-bonn.de
