function [rho_sigma,rho_epsilon,eta_sigma,sigma_bar]=par_transform_AR1(rho_sigma,rho_epsilon,eta_sigma,sigma_bar)
% [rho_sigma,rho_epsilon,eta_sigma,sigma_bar]=par_transform_AR1(rho_sigma,rho_epsilon,eta_sigma,sigma_bar)
% transforms parameters to unbounded domain
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

eta_sigma=log(eta_sigma);
rho_sigma= log((1-rho_sigma)./rho_sigma);
rho_epsilon= log((0.9999-rho_epsilon)./(rho_epsilon+0.9999)); %logistic transformation
sigma_bar=sigma_bar;