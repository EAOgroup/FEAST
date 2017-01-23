function [GasField] = GasFieldInitiator()

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

% This function creates a struct that defines a simulated gas field
% Inputs:
%   None
% Outputs:
%   GasField        Struct with the following fields:
%       LeakProductionRate  Rate leaks are created [leaks/component-day]
%       ComponentsPerSite   Number potential leak sources per well
%                           [components per well]
%       SiteCount           Number of wells to be simulated [wells]
%       MaxCount            Maximum number of wells allowed to be simulated
%                           by a single piece of technoloy. Represents the 
%                           number of wells in the field, rather than the 
%                           number of wells in the simulation. [wells]
%       SiteSpacing         Driving distance between wells [m]
%       WellArea            Area around each well that may have a leak
%                           [m^2]
%       H0max               Maximum height at which a leak may occur [m]
%       LeakData            Vector of possible leak sizes [g/s]
%       LeakSample          Vector of leak sizes in simulated field [g/s]
%       x0                  Vector of x positions of each leak
%       y0                  Vector of y positions of every leak [m]
%       H0                  Vector of heights of every leak [m]
%       Fonethird           Vector of buoyancy flux raised to the onethird
%                           for every leak [m^{4/3}/s]
%       LeakCount           Number of leaks [-]
%       ComponentCount      Number of Components [-]
%       NullRepairRate      Rate at which leaks are repaired through the
%                           null process [Repairs/leak-day]

%% User Inputs
%Number of leaks per well
LeaksPerWell = 6;

GasField.LeakProductionRate = 1e-5; %new leaks per component per day
%Number of valves and connectors per well (736659 total components/1138 wells in the Fort Worth study)
GasField.ComponentsPerSite=650;
%Number of wells to be simulated
GasField.SiteCount = 100; %1100;

%Maximum number of wells to be surveyed with a single capital investment
GasField.MaxCount = 6000;

%Driving distance between wells
GasField.SiteSpacing=700; %m

% Area overwhich a leak may be found per well. Based on satelite imagery.
GasField.WellArea = 100; %m^2

% Maximum leak height
GasField.H0max = 5; %m

%Number of leaking components at each well (arbitrary standard deviation)
LeaksInWell = round(2*randn(1,GasField.SiteCount+1)+LeaksPerWell);
LeaksInWell(LeaksInWell<0) = 0;

g=9.8; % g is the strength of gravity [m/s^2]
rhoAir = 1225; % density of air [g/m^3]
rhoMethane = 681; % density of methane at atmospheric pressure [g/m^3]

%% Leak Size Distribution
%load empirical leak distribution data (lbs/year methane). copied from the
%fort worth EmissionsData.xlsx file, sheet 'emissions data', column
%'methane', including components with no IR measurement.

dataStruct=importdata([pwd '/DataFiles/FortWorthScraped.mat']);
% Data from fort worth came from both traditional
% FID detection and IR camera detection. All of the components were
% surveyed with the IR camera, and 10% of the components were also surveyed
% with the FID. The FID found many leaks that the camera missed. To
% accurately represent the probability of finding a leak identified by the
% camera correctly, the leak sample only includes 10% of the leaks found by
% the camera.
GasField.LeakData = [dataStruct.FID; datasample(dataStruct.IR,round(length(dataStruct.IR)/10))]; %[g/s]
%data2=importdata([pwd '/LeakData/AllenSCFMMethaneScraped.mat')*1e5*2.8e-2*16/(8.314*273.15*60); %[g/s]

%Generate sets of leaks of various sizes. The list of leak sizes can be
%associated with specific wells using the LeaksInWell vector.
%Total Number of leaking components
LeakCount = sum(LeaksInWell);

GasField.LeakSample = datasample(GasField.LeakData,LeakCount);

%Now the leaks are distributed randomly about the center of the natural gas
%well.

%% Final Output
GasField.x0=sqrt(GasField.WellArea)*(rand(LeakCount,1)-0.5);
GasField.y0=sqrt(GasField.WellArea)*(rand(LeakCount,1)-0.5);
GasField.H0=GasField.H0max*rand(LeakCount,1);

Ffactor = g /pi*(1/rhoMethane-1/rhoAir); % Bouyancy flux [m^4/s^3]
GasField.Fonethird = (Ffactor * GasField.LeakSample).^(1/3);
ComponentCount = GasField.ComponentsPerSite * GasField.SiteCount ;
GasField.LeakCount = LeakCount; GasField.ComponentCount = ComponentCount;
GasField.NullRepairRate = GasField.LeakProductionRate*ComponentCount/LeakCount; % leaks repaired per leak per day
end
