function [ AtmosphereStruct] = AtmosphereInitiator(timeSteps)

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

% atmosphereInitiator creates a struct containing relevant atmospheric
% data.
% Inputs:
    % timeSteps     The number of time intervals to be simulted
% Outputs:
    % atmosphereStruct.windspeeds   a vector of randomnly generated windspeeds [m/s]
    
    
% Use ARPAE MONITOR wind data set
windData =importdata([pwd '/DataFiles/ARPAEWind.csv']);
windSpeedData = windData.data(:,2); %[m/s]
% Now use Fort Worth data for wind direction.
windData = importdata([pwd '/DataFiles/FortWorthWindData.csv']);
% windSpeedData = windData.data(:,1)*0.447; %[m/s]

% Generate set of windspeeds
windspeeds = datasample(windSpeedData,timeSteps);
AtmosphereStruct.u = windspeeds;

% Generate set of wind directions
windDirectionData = windData.data(:,3)*pi/180; %[radians]
AtmosphereStruct.theta = datasample(windDirectionData,timeSteps);

% Generate set of atmospheric stability classes
class = zeros(timeSteps,1);
class(windspeeds<2) = randi([1,2],size(windspeeds(windspeeds<2)));
class(windspeeds>=2 & windspeeds<3) = randi([1,3],size(windspeeds(windspeeds>=2 & windspeeds<3)));
class(windspeeds>=3 & windspeeds<5) = randi([2,4],size(windspeeds(windspeeds>=3 & windspeeds<5)));
class(windspeeds>=5) = randi([3,4],size(windspeeds(windspeeds>=5)));
AtmosphereStruct.class = class;
Ry = [0.23 0.18 0.11 0.08 0.058 0.04];
Rz = [0.18 0.12 0.07 0.05 0.035 0.025];

% Calculate Ry and Rz for use in the Gaussian plume model based on the
% stability class.
AtmosphereStruct.Ry = Ry(class); AtmosphereStruct.Rz = Rz(class);