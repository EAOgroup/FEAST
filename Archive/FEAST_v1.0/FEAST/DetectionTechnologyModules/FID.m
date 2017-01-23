function [FIDStruct, financeStruct, processStruct] = FID(GasField,timeStruct,varargin)

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

% This function creates structs that fully characterize a manually operated
% FID detection technique.
% Inputs:
%   GasField    Struct defining the natural gas field to be simulated
%   timeStruct  Struct containing the following elements:
%       time    Vector of times in the simulation
%       timeSteps Number of time steps in the simulation
%       endTime   The final time in the simuation
%       deltaT    Length of each timestep
%       varargin  Allows values to be input for sensitivity studies
% Outputs:
%   FIDStruct contains properties of the technology
%   financeStruct contains financial information about the FID
%   processStruct contains information about how gas fields will be surveyed.

%% varargin arguments

surveyTimeSensitivityFactor = 1;
surveyIntervalSensitivityFactor = 1;
lifetimeFactor = 1;

for i = 1:2:length(varargin)
    arg = varargin{i};
    if strcmpi(arg,'surveyTimeSensitivityFactor')
        surveyTimeSensitivityFactor = varargin{i+1};
    elseif strcmpi(varargin{i},'lifetimeSensitivityFactor')
        lifetimeFactor = varargin{i+1};
    elseif strcmpi(varargin{i},'surveyIntervalSensitivityFactor')
        surveyIntervalSensitivityFactor = varargin{i+1};
    end
end


%% FID Properties
FIDStruct.lifetime = 10*365*lifetimeFactor; %days

%% Process Details
processStruct.surveyInterval = 100*surveyIntervalSensitivityFactor; %days
processStruct.surveySpeed = 150; %components/person-hour
processStruct.driveSpeed = 15; %m/s
processStruct.setupTime = 0.1; %hours
processStruct.workTime = 5; % hours/day on average...5 hours per day on average corresponds to 35 hours per week.
processStruct.surveyTime = (GasField.ComponentsPerSite/processStruct.surveySpeed+GasField.SiteSpacing/processStruct.driveSpeed/3600+processStruct.setupTime)*surveyTimeSensitivityFactor; %hours
% timeFactor accounts for the finite simulation size. The effective capital cost is
% reduced in the simulation based on the ratio of the wells in the
% simulation to the number of wells that a single FID could survey.
processStruct.timeFactor = processStruct.surveyTime*GasField.SiteCount/(processStruct.surveyInterval*processStruct.workTime);

%% Financial Properties

Capital0 = 35000*processStruct.timeFactor; %dollars (covers 5k for FID and 30k for truck)
Maintenance0 = Capital0*0.1; %dollars/year

financeStruct.Capital = zeros(1,timeStruct.timeSteps); %dollars

lifeTimes=0;
while lifeTimes < timeStruct.endTime
    financeStruct.Capital(max(1,round(lifeTimes/timeStruct.deltaT)))=Capital0;
    lifeTimes = lifeTimes + FIDStruct.lifetime;
end

financeStruct.maintenance = ones(1,timeStruct.timeSteps)*Maintenance0*timeStruct.deltaT/365; % $

financeStruct.Labor = 100; %dollars/hour

% maintenance costs are estimated as 10% of capital per year
financeStruct.maintenance = ones(1,timeStruct.timeSteps)*Maintenance0*timeStruct.deltaT/365; % $

% surveyCost is the cost to survey all wells in the natural gas field
financeStruct.surveyCost = financeStruct.Labor*GasField.SiteCount*processStruct.surveyTime;

% Find cost is the cost of searching for leaks
financeStruct.findCost = zeros(1,timeStruct.timeSteps);
financeStruct.findCost(mod(timeStruct.time,processStruct.surveyInterval)<timeStruct.deltaT) = financeStruct.surveyCost;