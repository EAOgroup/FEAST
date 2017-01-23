function [ leak ] = LeakCreater( Ncomponents, deltaT, LeakProductionRate, LeakData )

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

%Creates new leaks in a gas field
%   Inputs
%       Ncomponents     Number of components that may be leaking
%       deltaT          Period of time to be simulated [days]
%       LeakProductionRate Rate at which new leaks are formed [leaks per
%                           component per day]
%       LeakData        Distribution of leak sizes to draw from [g/s]
%   Outputs:
%       leak            struct with the following fields:
%           Q               vector of sizes of every leak [g/s]
%           x               vector of x positions of every leak [m]
%           y               vector of y positions of every leak [m]
%           H               vector of heights of every leak [m]
%           Fonethird       vector of buoyancy flux raised to the onethird
%                           for every leak [m^{4/3}/s]
%           leaksDetected   vector initiated to zeros. Intended to identify
%                           leaks that have been detected but not repaired.

% Calculate constants
g=9.8; % g is the strength of gravity [m/s^2]
rhoAir = 1225; % density of air [g/m^3]
rhoMethane = 681; % density of methane at atmospheric pressure [g/m^3]
Ffactor = g /pi*(1/rhoMethane-1/rhoAir); % Bouyancy flux [m^4/s^3]

%Randomly choose the total number of leaks
List = rand(1,Ncomponents);
LeakCount = sum(List<LeakProductionRate * deltaT);

% Randomly pull leak sizes from the input data
leak.Q = datasample(LeakData,LeakCount);

%Randomly define other leak properties
leak.x=10*(rand(LeakCount,1)-0.5);
leak.y=10*(rand(LeakCount,1)-0.5);
leak.H=5*rand(LeakCount,1);
leak.Fonethird = (Ffactor * leak.Q).^(1/3);
leak.leaksDetected = zeros(LeakCount,1);
end

