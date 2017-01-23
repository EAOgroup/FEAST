function [camQmin] = DetectionLimit(camStruct,gasField,varargin)

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

% detectionLimit determines the smallest leak that can be detected by an IR
% camera based technology, based on parameters of the camera and optimal
% atmospheric conditions.

%% Inputs:
%   camStruct       Structure containing at least the following camera properties:
%       R0              1x3 Position of camera [m]
%       n               1x3 direction that the camera points
%       top             1x3 direction of the top of the camera
%       right           1x3 direction of the right side of the camera
%       npoints         square root of the number of pixels in the
%                       camera
%       minDetection    minimum detectable concentration pathlength [g-m/m^3]
%       FOV1            Field of view in the left-right orientation [radians]
%       FOV2            Field of view in the top-bottom orientation [radians]
%       theta           Array of polar angles associated with each pixel
%                       [radians]
%       phi             Array of azimuth angles associated with each pixel
%                       [radians]
%       minPixelCount   minimum fraction of pixels with a signal above minDetection
%                       for the leak to be found.
%   GasField        Structure containing at least the following fields:
%       H0max           Maximum allowable height for a leak [m]
%       LeakData        Leak distribution [g/s]
%   varargin        Optional input to state how camera height is defined

%% Read optional input
if length(varargin) == 2
    if ~strcmpi(varargin{1},'heightdefinition')
        display('invalid property name in detectionLimit');
        return
    end
    if ~strcmpi(varargin{2},'leak') && ~strcmp(varargin{2},'ground')
        display('invalid height definition in LeakDetector');
        return
    end
    HeightDef = varargin{2};
else
    HeightDef = 'ground';
end

%% Perorm calculation
% Use ARPAE MONITOR wind data set
winddata =importdata([pwd '/DataFiles/ARPAEWind.csv']);
windspeeddata = winddata.data(:,2); %[m/s]
atmosphere.u = min(windspeeddata);
atmosphere.class = 6;

leak.H = gasField.H0max;
leak.x = 0;
leak.y = 0;
g=9.8; % g is the strength of gravity [m/s^2]
rhoAir = 1225; % density of air [g/m^3]
rhoMethane = 681; % density of methane at atmospheric pressure [g/m^3]
Ffactor = g /pi*(1/rhoMethane-1/rhoAir); % Bouyancy flux [m^4/s^3]

leak.Q = min(gasField.LeakData);
leak.Fonethird = (Ffactor * leak.Q).^(1/3);
iterations = 0;
while LeakDetector(camStruct,leak,atmosphere,'HeightDefinition',HeightDef) ~= 1
    leak.Q = 2*leak.Q;
    leak.Fonethird = (Ffactor * leak.Q).^(1/3);
end

camQmin = leak.Q/2;