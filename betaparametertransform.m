function [a,b]=betaparametertransform(betamean,betavariance)
% function [a,b]=betaparametertransform(betamean,betavariance)
% transforms mean and variance of beta distribution to a and b used in
% Matlab
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

b=(betamean-2*betamean^2+betamean^3-betavariance+betamean*betavariance)/betavariance;
a=-betamean*b/(betamean-1);