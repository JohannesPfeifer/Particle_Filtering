function [logpost,x_hat,x_std,Eff_particle]=PF_caller(draw,observable_series,num_sim_filter,num_sim_smoother,smoother_dummy,name)
% [loglik, flag, xhat_PF0,x_std_PF0,Eff_particle]=Particle_smoother(draw,observable_series,num_sim,graph_dummy,name)
%Inputs:
%   draw                [4 by 1] vector                 vector of estimated parameters
%   observable_series   [T by 1] vector                 observed data
%   num_sim_filter      scalar                          number of particles
%   num_sim_smoother    scalar                          number of particles for smoother
%   smoother_dummy      scalar                          dummy to call smoother
%   name                string                          name of the file under which to save
%                                                       from smoother
%
% Output
%  logpost              scalar            The log-posterior
%  x_hat                [T by 1] vector   the posterior state estimate
%                                           x_t|T
%  x_std                [T by 1] vector   standard deviation of the posterior
%  Eff_particle     	[T by 1] vector   A measure of the effective sample size
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


if nargin<6
   name='smootherresults';
end

%initialize variables
T=length(observable_series);
x_hat=NaN(T,1);
x_std=NaN(T,1);
Eff_particle=NaN(T,1);

rho_sigma=draw(1);
rho=draw(2);
eta_sigma=draw(3);
sigma_bar=draw(4);

[rho_sigma,rho,eta_sigma,sigma_bar]=par_retransform_AR1(rho_sigma,rho,eta_sigma,sigma_bar);

[priorval,alarm]=evaluate_prior_AR1([rho_sigma,rho,eta_sigma,sigma_bar]);
if alarm
    logpost=Inf;
    return
end

% *********************** Filtering/Smoothing ***************************************

x0 = sigma_bar+eta_sigma*randn(num_sim_filter,1); % The starting values for filtering from initial distribution

[LogL,Eff_particle,x_hat,x_std] = PF_basefunction_with_smoother_AR1(@PF_AR1_density,x0,observable_series,num_sim_filter,num_sim_smoother,rho_sigma,rho,sigma_bar,eta_sigma,name,smoother_dummy);
% add prior
logpost=sum(LogL)+priorval;
