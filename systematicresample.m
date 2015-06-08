function N_index= systematicresample(weights,randnr,num_sim)
% This function implements systematic resampling, based on Arulampalam/Maskell/Gordon/Clapp (2002)
%   "A Tutorial on Particle Filters for Online Nonlinear/Non-Gaussian
%   Bayesian Tracking", IEEE Transactions on Signal Processing, 50(2)  
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


% Calculating the cdf for the weights
cdf_weights = cumsum(weights);

% A starting value for u(1)
constant  = 1/num_sim;
u=(randnr*constant:constant:1-(1-randnr)*constant)';
i         = 1;
N_index    = zeros(num_sim,1);
indexcount = 1;
for j=1:num_sim
   % Move along the cdf
   while u(j) > cdf_weights(i)
      if (i < num_sim) 
          i = i + 1;
      end
      indexcount=indexcount+1;
   end
  N_index(j,1)=indexcount;
end
