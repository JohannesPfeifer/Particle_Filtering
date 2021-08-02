function [priorval,alarm]=evaluate_prior_AR1(par)
%compute the prior distribution for
%par=[rho_sigma,rho,eta_sigma,sigma_bar] of a generic AR1 stochastic
%volatility process
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

alarm=0;
priorval=0; % intialize prior

% rho_sigma~b(0.9,0.1^2)
[a b]=betaparametertransform(0.9,0.1^2);
ii=[1];
prioradd=log(betapdf(0.999^(-1)*par(ii)',a,b));
if ~isinf(prioradd)
  priorval=priorval+sum(prioradd);
else
  alarm=1;  
  return
end

% rho ~ U(-0.9999,0.9999 )
ii=[2];
if par(ii)>=-0.9999 && par(ii)<=0.9999 
  priorval=priorval-log(1.9998);
else
  alarm=1;  
  return
end

%eta_sigma ~ Gamma(0.5,0.1^2);
[a, b]=gammaparametertransform(0.5,0.1^2);
ii=[3];
prioradd=sum(log(gampdf(par(ii)',a,b)));
if ~isinf(prioradd)
  priorval=priorval+prioradd;
else
  alarm=1;  
  return
end

%sigma_bar ~ U(-11,-3)
ii=[4];
if par(ii)>=(-11) && par(ii)<-3
  priorval=priorval-log(8);
else
  alarm=1;  
  return
end

end
