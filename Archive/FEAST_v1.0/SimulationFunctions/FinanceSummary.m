function [ Capital, Maintenance, Finding ] = FinanceSummary(financeStruct)

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

%This function takes results from the input arguments and combines them into arrays 
%   Inputs:
%       financeStruct must be an array of structs.The structs
%       must have the following vector fields of the same length:
%            findCost
%            maintenance
%            Capital
%   Outputs:
%       Capital      Array of capital costs through time
%       Maintenance  Array of maintenance costs through time
%       Finding      Array of finding costs through time

% Example of use:
% [ Capital Maintenance Finding ] = financeSummary(financeStruct)

Ntech = length(financeStruct);

Ntimes = length(financeStruct{1}.Capital);
Capital = zeros(Ntech,Ntimes); Maintenance = zeros(Ntech,Ntimes); Finding = zeros(Ntech,Ntimes);

for i=1:Ntech
    Capital(i,:) = financeStruct{i}.Capital;
    Maintenance(i,:) = financeStruct{i}.maintenance;
    Finding(i,:) = financeStruct{i}.findCost;
end

