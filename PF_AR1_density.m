function [density]= PF_AR1_density(y,x_t,rho_epsilon1,num_sim,dimy,t)
% function [density]= PF_AR1_density(y,x_t,rho_epsilon1,num_sim,dimy,t)
% This function evaluates the density p(y(t)|x(t)) used in a Standard Particle filter
%
%Inputs:
%   y                   vector      observed data
%   x                   vector      state x_t
%   rho_epsilon1        scalar      AR1-coefficient level equation
%   num_sim             scalar      Number of particle used
%   t                   scalar      current observation used for reading
%                                   from y
%
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


density = zeros(num_sim,1); 
cons_density  = (2*pi)^(-dimy*0.5)./exp(x_t);  

y_t = y(t)-rho_epsilon1*y(t-1);
w = bsxfun(@rdivide,y_t,exp(x_t));  

density(:,1) = exp(-0.5*w.^2);       
density = cons_density.*density;
