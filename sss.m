function [y_sss,y] = sss(dr, long, iorder)

% function y = irf(dr, long, iorder)
% Computes the stochastic steady state
% 
% INPUTS
%    dr:     structure of decisions rules for stochastic simulations
%    long:   number of periods of simulation
%    iorder: first or second order approximation
%
% OUTPUTS
%    y_sss:      vector of stochastic steady states
%        
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2010 Dynare Team
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

global M_ oo_ options_

if M_.dynare_version == '4.5.7'
    if M_.maximum_lag >= 1
    temps = repmat(dr.ys,1,M_.maximum_lag);
else
    temps = zeros(M_.endo_nbr, 1); % Dummy values for purely forward models
end

ex1     = zeros(long,M_.exo_nbr);
y       = simult_(temps,dr,ex1,iorder);
y_sss   = y(:,end);
    
else

if M_.maximum_lag >= 1
    temps = repmat(dr.ys,1,M_.maximum_lag);
else
    temps = zeros(M_.endo_nbr, 1); % Dummy values for purely forward models
end

ex1     = zeros(long,M_.exo_nbr);
y       = simult_(M_, options_, temps,dr,ex1,iorder);
y_sss   = y(:,end);
end
end
