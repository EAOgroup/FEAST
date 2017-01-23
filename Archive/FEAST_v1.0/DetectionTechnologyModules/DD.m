function [snifferStruct, financeStruct, processStruct] = DD(GasField,timeStruct,varargin)

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

% This model simulates the performance of a network of distributed
% sniffers.
% Inputs:
%   GasField    Struct defining the natural gas field to be simulated
%   timeStruct  Struct containing the following elements:
%       time    Vector of times in the simulation
%       timeSteps Number of time steps in the simulation
%       endTime   The final time in the simuation
%       deltaT    Length of each timestep
% Outputs:
%   snifferStruct Contains variables defining the technology
%   financeStruct contains financial information about the sniffer
%   processStruct contains information about how gas fields will be surveyed

%% varargin arguments
% surveyTimeFactor
% lifetimeFactor

surveyTimeSensitivityFactor = 1;
lifetimeFactor = 1;

for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'surveyTimeSensitivityFactor')
        surveyTimeSensitivityFactor = varargin{i+1};
    elseif strcmpi(varargin{i},'lifetimeSensitivityFactor')
        lifetimeFactor = varargin{i+1};
    end
end

%% Detector properties
%Sensitivity is the trigger sensitivity for the sniffer
snifferStruct.Sensitivity = .01; %g/m^3 (multiply by 1500 to get ppm)
snifferStruct.x = [-3.5 0 3.5 0 ];
snifferStruct.y = [-10, -10, -10, 10 ];
snifferStruct.z = [2, 3, 4, 4];
snifferStruct.N = length(snifferStruct.x); % sniffers per well
snifferStruct.lifetime = 5*365*lifetimeFactor; %days

%% Process details
processStruct.repairInterval = 50; % days
processStruct.surveySpeed = 500; %component/person-hour
processStruct.driveSpeed = 15; %m/s
processStruct.setupTime = 0.5; %hours
processStruct.locationTime = (GasField.ComponentsPerSite/processStruct.surveySpeed/snifferStruct.N+GasField.SiteSpacing/processStruct.driveSpeed/3600+processStruct.setupTime)/snifferStruct.N*surveyTimeSensitivityFactor; % hours to locate one leak
processStruct.workTime = 5; % hours/day on average...5 hours per day on average corresponds to 35 hours per week.
processStruct.timeFactor = processStruct.locationTime*GasField.SiteCount/(processStruct.repairInterval*processStruct.workTime);

%% Financials
SnifferCost = 500; %$ per sniffer
financeStruct.SnifferCapital = SnifferCost*snifferStruct.N; % $ per well
financeStruct.CameraCapital = 90000;%$
financeStruct.TruckCapital = 30000; %$
Capital0 = financeStruct.SnifferCapital*GasField.SiteCount+processStruct.timeFactor*(financeStruct.CameraCapital + financeStruct.TruckCapital); %$
Maintenance0 = Capital0*0.1;
financeStruct.maintenance = ones(1,timeStruct.timeSteps)*Maintenance0*timeStruct.deltaT/365;
financeStruct.Labor = 100; %$/hour
financeStruct.locationCost = financeStruct.Labor*processStruct.locationTime; %$/leak

% Calculate capital costs vector
financeStruct.Capital = zeros(1,timeStruct.timeSteps); %dollars

lifeTimes=0;
while lifeTimes < timeStruct.endTime
financeStruct.Capital(max(1,round(lifeTimes/timeStruct.deltaT)))=Capital0;
lifeTimes = lifeTimes + snifferStruct.lifetime;
end

% The finding costs must be calculated as leaks are found. The vector is
% initialized here.
financeStruct.findCost = zeros(1,timeStruct.timeSteps);