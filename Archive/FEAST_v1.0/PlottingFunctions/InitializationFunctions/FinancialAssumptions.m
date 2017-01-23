function [financeStruct] = FinancialAssumptions()

%-------------------------------------------------------------------------------------------------------------------
% This file is part of Fugitive Emissions Abatement Simulation Toolkit aka FEAST
% Copyright {2016} {Chandler E. Kemp; Arvind P. Ravikumar; Adam R. Brandt} 

% FEAST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% FEAST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------------------------------------------------

%   FinancialAssumptions creates a struct that stores financial  data to be
%   used in the simulation.

% Price of gas
financeStruct.gasPrice = 2e-4; %dollars/g (2e-4 dollars/g = $5/mcf-methane at STP-273.15 K, 1 bar)
% Discount rate for NPV analysis.
financeStruct.discountRate = 0.08;