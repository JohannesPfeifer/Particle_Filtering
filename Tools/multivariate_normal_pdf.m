function density = multivariate_normal_pdf(X,Mean,Sigma_upper_chol,n)
% Evaluates the density of a multivariate gaussian, with expectation Mean
% and variance Sigma_upper_chol'*Sigma_upper_chol, at X.
%
%
% INPUTS
%
%    X                  [double]    1*n vector
%    Mean               [double]    1*n vector, expectation of the multivariate random variable.
%    Sigma_upper_chol   [double]    n*n matrix, upper triangular Cholesky decomposition of Sigma (the covariance matrix).
%    n                  [integer]   dimension.
%
% OUTPUTS
%    density            [double]    density
%
% SPECIAL REQUIREMENTS

% Copyright (C) 2003-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
density = (2*pi)^(-.5*n) * ...
          prod(diag(Sigma_upper_chol))^(-1) * ...
          exp(-.5*(X-Mean)*(Sigma_upper_chol\(transpose(Sigma_upper_chol)\transpose(X-Mean))));