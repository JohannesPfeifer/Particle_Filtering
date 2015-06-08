function [rho_sigma,rho,eta_sigma,sigma_bar]=par_retransform_AR1(rho_sigma,rho,eta_sigma,sigma_bar)
% [rho_sigma,rho,eta_sigma,sigma_bar]=par_retransform_AR1(rho_sigma,rho,eta_sigma,sigma_bar)
% retransform parameter from unbounded support to original values
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

if nargin==1
    if size(rho_sigma,2)==4
       sigma_bar=rho_sigma(:,4);
       rho=rho_sigma(:,2);
       eta_sigma =rho_sigma(:,3);
       rho_sigma=rho_sigma(:,1);
    else
        error('Wrong dimensions')
    end
end

eta_sigma=exp(eta_sigma);
rho_sigma= 1./(1+exp(rho_sigma));
rho=  (0.9999+exp(rho)*(-0.9999))./(1+exp(rho));
sigma_bar=sigma_bar;

if nargout==1
    rho_sigma=[rho_sigma,rho,eta_sigma,sigma_bar];
end

