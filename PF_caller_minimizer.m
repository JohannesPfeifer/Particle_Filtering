function [minus_logpost,x_hat,x_std,Eff_particle]=PF_caller_minimizer(varargin)
% [minus_logpost, flag, xhat_PF0,x_std_PF0,Eff_particle]=PF_caller_minimizer(varargin)
%Inputs:
%   varargin            cell array        all required input arguments
%
% Output
%  logpost              scalar            The log-posterior
%  x_hat                [T by 1] vector   the posterior state estimate
%                                           x_t|T
%  x_std                [T by 1] vector   standard deviation of the posterior
%  Eff_particle     	[T by 1] vector   A measure of the effective sample size
%
% Copyright (C) 2013-2015 Benjamin Born + Johannes Pfeifer
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


[logpost,x_hat,x_std,Eff_particle]=PF_caller(varargin{:});
minus_logpost=-logpost;