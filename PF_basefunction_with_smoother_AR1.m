function [LogL,Eff_particle,x_hat,x_std] = PF_basefunction_with_smoother_AR1(PF_density_computation_function,...
    x0,y,num_sim_filter,num_sim_smoother,rho_sigma,rho_1,sigma_bar,eta_sigma,name,smoother_dummy)
%   This function calculates the standard Particle Filter for the system
%   y(t)   =rho*y(t)+ e^sigma(t)*u(t)
%   sigma(t+1) =(1-rho_sigma)*sigma_bar+rho_sigma*sigma(t) + eta_sigma*eps(t+1) 
% 
% Implemented Filter according to Arulampalam/Maskell/Gordon/Clapp (2002)
%   "A Tutorial on Particle Filters for Online Nonlinear/Non-Gaussian
%   Bayesian Tracking", IEEE Transactions on Signal Processing, 50(2)  
% 
% Implemented smoother according to Godsill/Doucet/West(2004) "Monte Carlo smoothing for nonlinear time series", 
%   Journal of the American Statistical Association, 2004, 99, 156-168

% Input
%   PF_density_computation_function - function evaluating the importance sampling weights
%   observable_series   
%   x0                  [num_sim_filter by 1] vector    Initial state vector of dimension
%   y                   [T by 1] vector                 observed data
%   num_sim_filter      scalar                          number of particles
%   num_sim_smoother    scalar                          number of particles for smoother
%   rho_sigma,rho_1,sigma_bar,eta_sigma scalars         Parameters of the LOM estimated with filter.
%   smoother_dummy      scalar                          dummy to call smoother
%   name                string                          name of the file under which to save
% 
% Output
%  LogL                 [T by 1] vector   The contribution to the log-likelihood function 
%                                         for each observation
%  Eff_particle     	[T by 1] vector   A measure of the effective sample size
%  x_hat                [T by 1] vector   the posterior state estimate
%                                           x_t|T
%  x_std                [T by 1] vector   standard deviation of the posterior

%
% Copyright (C) 2013-2014 Benjamin Born + Johannes Pfeifer
%
% This is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% See <http://www.gnu.org/licenses/>.


% >>>>>>>>>>>>>>>>>>>>>>>>>>> INITIALIZATIONS <<<<<<<<<<<<<<<<<<<<<<<<<<
dimy=size(y,2);
dimx=size(x0,2);

T_obs         = size(y,1);      % The sample size
%draw random numbers used
shocks = randn(num_sim_filter,T_obs);
randnr = rand(T_obs,1);

% Setting dimensions
x_hat           = zeros(T_obs,dimx);   % Matrix for the posterior state estimates
x_std           = zeros(T_obs,dimx);   % Matrix for the posterior state estimate standard deviation
if smoother_dummy
    x_std          = zeros(T_obs,dimx);   % Matrix for the posterior state standard deviation
end
LogL           = zeros(T_obs,1);      % Vector for storing the contributions to the log-likelihood fct
Eff_particle   = zeros(T_obs,1);      % Vector for the effective sample size
x_resampled_store = zeros(dimx,num_sim_filter,size(y,2));
x_simulated_store = zeros(dimx,num_sim_filter,size(y,2));
weights_store= zeros(dimx,num_sim_filter,size(y,2));
weights_store(1,:,1)= 1/num_sim_filter;

% The first observation 
x_resampled                = x0;
x_resampled_store(1,:,1)   = x0;
x_simulated_store(1,:,1)   = x0;
x_hat(1,1)                 = mean(x0);


% >>>>>>>>>>>>>>>>>>>>>>>>>>>> FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<
for t=2:T_obs
   % We draw from the state transition distribution using draws from the 
   % previous state distribution (x_resampled) and the shocks (shocks(:,t))
   x_simulated = PF_stoch_vol_LOM(x_resampled,rho_sigma,sigma_bar,eta_sigma,shocks(:,t));

   % We evaluate the importance sampling weights which are given by p(y|x) which we call "density" 
   [density]   = feval(PF_density_computation_function,y,x_simulated,rho_1,num_sim_filter,dimy,t);

   % The sum of density
   sum_density = sum(density);
   
   % We compute the normalized importance sampling weights
   if sum_density < 1e-12
       LogL=-Inf;
       disp(t)
       %error('No effective particles')
       return;
   else
      weights      = density/sum_density;
   end
   weights_store(1,:,t)    = weights;
   % Contribution to the log-likelihood function 
   LogL(t,1) = log(sum_density/num_sim_filter);

   % The effective sample size
   Eff_particle(t,1) = 1/sum(weights.^2);
   
   % disp(Eff_particle(t,1))
   % We resample the swarm of particles using systematic resampling.

   N_indices  = systematicresample(weights,randnr(t,1),num_sim_filter); 
   x_resampled=x_simulated(N_indices);
   x_resampled_store(1,:,t)=x_resampled;
   x_simulated_store(1,:,t)=x_simulated;

   % Saving the posterior state estimate
   x_hat(t,:) = mean(x_resampled,1);
   x_std(t,:) = std(x_resampled,0,1);
end

if smoother_dummy
    %initialize random numbers for drawing states
    rand_numbers=rand(num_sim_smoother,T_obs);

    Eff_particle_smoother=zeros(num_sim_smoother,T_obs);
    x_smoother_store=zeros(num_sim_smoother,T_obs);

    for traject_iter=1:num_sim_smoother;
        %choose last state
        cum_prob=cumsum(weights_store(1,:,T_obs));
        x_tilde_plus1=x_simulated_store(1,sum(cum_prob<rand_numbers(traject_iter,T_obs))+1,T_obs);
        x_smoother_store(traject_iter,T_obs)=x_tilde_plus1;
        %% smoother loop
        for t=T_obs-1:-1:1
            % draw from time T distribution
           x_filter=x_simulated_store(1,:,t); 

           raw_density  = Smoother_state_density(x_tilde_plus1,x_filter,rho_sigma,sigma_bar,eta_sigma,num_sim_filter,dimx);

           density= raw_density.*weights_store(1,:,t)';
           % The sum of density
           sum_density = sum(density);

           % We compute the normalized importance sampling weights
           if sum_density < 1e-12
               disp(t)
               error('Zero density')
           else
              weights      = density/sum_density;
           end

           % The effective sample size
           Eff_particle_smoother(traject_iter,t) = 1/sum(weights.^2);

            cum_prob=cumsum(weights);
           x_t_smoothed=x_filter(1,sum(cum_prob<rand_numbers(traject_iter,t))+1);

           % Saving the posterior state estimate
           x_smoother_store(traject_iter,t)=x_t_smoothed;
           %set for next period
           x_tilde_plus1=x_t_smoothed;
        end
    end
    xhat_smoother = median(x_smoother_store,1);
    x_std_smoother = std(x_smoother_store,0,1);
    x_percentile_smoother = quantile(x_smoother_store,[0.025 0.975]);

    vola_resids=zeros(T_obs,num_sim_smoother);
    vola_resids(2:T_obs,:)=x_smoother_store(:,2:T_obs)'-(1-rho_sigma)*sigma_bar-rho_sigma*x_smoother_store(:,1:end-1)';
    level_resids=zeros(T_obs,num_sim_smoother);
    level_resids(2:T_obs,:)=bsxfun(@rdivide,(y(2:T_obs)-rho_1*y(1:T_obs-1))*ones(1,num_sim_smoother),exp(x_smoother_store(:,2:end)'));

    save(name,'y','rho_sigma','rho_1','sigma_bar','eta_sigma','x_percentile_smoother','x_std_smoother','xhat_smoother','Eff_particle_smoother','vola_resids','level_resids','x_smoother_store','x_resampled_store'); 
end