function density = Smoother_state_density(x_tilde_plus1,x,rho_sigma,sigma_bar,sigeta,num_sim,dimx)
% This function evaluates the state density p(x(t+1)|x(t)) used in the
% Smoother
%
%Inputs:
%   x_tilde_plus1               vector      x(t+1) chosen in previous period
%   x                           vector      state x_t from filtering
%   rho_sigma,sigma_bar,sigeta  scalars     coefficient for state transition equation 
%   num_sim                     scalar      Number of particle used
%   dimx                        scalar      number of states
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
cons_density  = (2*pi)^(-dimx*0.5)./sigeta;  % only valid of dimx=1

resid=x_tilde_plus1-(1-rho_sigma)*sigma_bar-rho_sigma*x;

w = bsxfun(@rdivide,resid,sigeta);  
   % The unscaled value of density
density(:,1) = exp(-0.5*w.^2);       
% Scaling the density by the constant
density = cons_density.*density;