function x_tplus1 = PF_stoch_vol_LOM(x_t,rho_sigma,sigma_bar,sigeta,shocks)
%function x_cup = PF_stoch_vol_LOM(x_t,rho_sigma,sigma_bar,sigeta,shocks)
% This function draws x_tplus1 (i.e. x(t+1)) from the state transition
% function given x_t and shocks
%
% Inputs:
%     - x_t   vector of states at time t
%     - rho_sigma,sigma_bar,sigeta    parameters used in LOM
%     - shocks    vector of stochastic shocks for transition
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

x_tplus1 = (1-rho_sigma)*sigma_bar+rho_sigma*x_t+sigeta*shocks;

