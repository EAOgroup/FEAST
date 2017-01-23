function [NullStruct, financeStruct, processStruct] = Null(GasField,timeStruct,varargin)

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

% This function creates to structs that fully characterize a manually operated
% IR detection technique.
% Inputs:
%   GasField    Struct defining the natural gas field to be simulated
%   timeStruct  Struct containing the following elements:
%       time    Vector of times in the simulation
%       timeSteps Number of time steps in the simulation
%       endTime   The final time in the simuation
%       deltaT    Length of each timestep
%       varargin  Allows values to be input for sensitivity studies

% Outputs:
%   NullStruct contains properties of the Null repair program
%   financeStruct contains financial information about the camera
%   processStruct contains information about how gas fields will be surveyed.

%% Null Finding Properties
% The null repair rate defaults to the same setting as the gas field, but
% can be set to a different value.
NullStruct.RepairRate = GasField.NullRepairRate; %GasField.LeakProductionRate*GasField.ComponentsPerSite * GasField.SiteCount / length(GasField.LeakSample);
%% Process details
% No process details to define here.
processStruct.field = [];
%% Financial Properties
% Capital costs are zero in this module.
financeStruct.Capital = zeros(1,timeStruct.timeSteps); %dollars

% maintenance costs are zero in the null module
financeStruct.maintenance = zeros(1,timeStruct.timeSteps); % $

% Find cost is zero in the null module
financeStruct.findCost = zeros(1,timeStruct.timeSteps);
