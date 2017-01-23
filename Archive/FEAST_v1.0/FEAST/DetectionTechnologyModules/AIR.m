function [camStruct, financeStruct, processStruct] = AIR(GasField,timeStruct,varargin)

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

% This function creates to structs that fully characterize a drone based
% automated IR detection technique.
% Inputs:
%   GasField    Struct defining the natural gas field to be simulated
%   timeStruct  Struct containing the following elements:
%       time    Vector of times in the simulation
%       timeSteps Number of time steps in the simulation
%       endTime   The final time in the simuation
%       deltaT    Length of each timestep
%       varargin  Allows values to be input for sensitivity studies

% Outputs:
%   camStruct contains properties of the IR camera
%   financeStruct contains financial information about the camera
%   processStruct contains information about how gas fields will be surveyed.

%% varargin arguments

surveyTimeSensitivityFactor = 1;
lifetimeFactor = 1;

for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'surveyTimeSensitivityFactor')
        surveyTimeSensitivityFactor = varargin{i+1};
    elseif strcmpi(varargin{i},'lifetimeSensitivityFactor')
        lifetimeFactor = varargin{i+1};
    end
end

%% Camera Properties
%%%%%%%%%%%%%%%%%%%%%Camera Parameters%%%%%%%%%%%%%%%%%%%%
FOV1=5.5*pi/180;                     %Horizontal field of view angle [radians]
FOV2=4.4*pi/180;                     %Vertical field of view angle [radians]
npoints=100;                         %pixels in the model in x and y directions.
minDetection=0.4;                    %minimum detectable concentration pathlength (g-m/m^3)
minPixelCount = 0.1;                 %minimum fraction of pixels with a signal above minDetection
%for the leak to be found.

%%%%%%%%%%%%%%%%%%%%%%%Note on Coordinates%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The coordinates here are defined such that x increases in the direction of the wind,
% z is equal to the distance above the ground, and y completes a right hand coordinate system.

%%%%%%%%%%%%%%Camera Location and orientation%%%%%%%%%%%%%%
camHeight = 20;                      %Camera height [m]
R0=[0,0,camHeight];                  %Camera position [m]
n=[0,0,-1];                          %Direction that the camera points [no units]
top = [0,1,0];                       %Direction that the top of the camera faces [no units]

camStruct = CreateCameraStruct('FOV1',FOV1,'FOV2',FOV2,'npoints',npoints,'minDetection',minDetection,'R0',R0,'n',n,'top',top,'minPixelCount',minPixelCount);
camStruct.lifetime = 3*365*lifetimeFactor; %days
camStruct.Qmin = DetectionLimit(camStruct,GasField,'HeightDefinition','ground');
%% Process details
processStruct.surveyInterval = 14; %days
processStruct.drivingSpeed = 15; %m/s
processStruct.surveySpeed = 5;%m/s
processStruct.setupFactor = 1.3; % allocates 20% of flight time to landing, refueling, and launching quadcopter
processStruct.workTime = 5; % hours/day on average...5 hours per day on average corresponds to 35 hours per week.
processStruct.surveyTime = ((GasField.WellArea/(camStruct.imageWidth*processStruct.surveySpeed/2)+GasField.SiteSpacing/processStruct.drivingSpeed)*processStruct.setupFactor/3600)*surveyTimeSensitivityFactor; % hours/well
% timeFactor accounts for the finite simulation size. The effective capital cost is
% reduced in the simulation based on the ratio of the wells in the
% simulation to the number of wells that a single camera could survey.
processStruct.timeFactor = max(processStruct.surveyTime*GasField.SiteCount/(processStruct.surveyInterval*processStruct.workTime),GasField.SiteCount/GasField.MaxCount);
%% Financial Properties
% Based on $100,000 for two A6700SC cameras, $30k truck, $50k UAV (based on
% quote for 11.5 kg payload + fuel "Penguin B UAV platform"), 10k filters
% and 2500 for computing power.
Capital0 = 193000*processStruct.timeFactor; %dollars
Maintenance0 = 0.1 * Capital0; %dollars/year

financeStruct.Capital = zeros(1,timeStruct.timeSteps); %dollars

lifeTimes=0;
while lifeTimes < timeStruct.endTime
    financeStruct.Capital(max(1,round(lifeTimes/timeStruct.deltaT)))=Capital0;
    lifeTimes = lifeTimes + camStruct.lifetime;
end

%LaborRate is the cost of labor
financeStruct.Labor = 100; %dollars/hour

% maintenance costs are estimated as 10% of capital per year
financeStruct.maintenance = ones(1,timeStruct.timeSteps)*Maintenance0*timeStruct.deltaT/365; % $

% surveyCost is the cost to survey all wells in the natural gas field
financeStruct.surveyCost = financeStruct.Labor*GasField.SiteCount*processStruct.surveyTime;

% Find cost is the cost of searching for leaks
financeStruct.findCost = zeros(1,timeStruct.timeSteps);
financeStruct.findCost(mod(timeStruct.time,processStruct.surveyInterval)<timeStruct.deltaT) = financeStruct.surveyCost;