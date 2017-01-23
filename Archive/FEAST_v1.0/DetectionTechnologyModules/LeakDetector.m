function [ detected ] = LeakDetector( cam, leak, atmosphere,varargin)

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

% This function determines which leaks are detected using an IR camera
% Inputs:
%   cam     A struct with the following fields:
%       R0     3D position of the camera [m]
%       n      direction that the camera points
%       top    direction that the top of the camera points
%       CamRight    direction that the right side of the camera points
%       npoints     square root of the number of pixels to be simulated
%       minDetection    minimum detectable concentration pathlength [g-m/m^3]
%       theta   array of polar angles associated with each pixel in the camera
%               [radians]
%       phi    array of azimuth angles associated with each pixel in the
%               camera [radians]
%       minPixelCount   the minimum number of pixels that must be have a signal
%        above cam.minDetection in order for the plume to be detected. Must be
%        between 0 and 1.
%   leak    A struct with the following fields:
%       Q      vector containing the flux of each leak [g/s]
%       H     height of the leak [m],
%       x, y define the location of each leak [m]
%       Fonethird   Buoyancy flux raised to the one-third power [(m^4/s^3)^1/3]
%   atmosphere    A struct with the following fields:
%       u      windspeed [m/s]
%       class  Defines the stability class of the atmosphere surrounding the
%        leak.
%       varargin    defines whether height is defined relative to the leak or to
%                   the ground
% Outputs: 
%   detected    A vector of ones and zeros showing which leaks were
%               detected.
%% Pull commonly used variables from structs
Q = leak.Q; H0 = leak.H; x0 = leak.x; y0 = leak.y; Fonethird = leak.Fonethird;
u = atmosphere.u; class = atmosphere.class;

%% Initiate variables
PixelCount = 0;
PixelCountOld = 0;
detected = 0;
minDetection = 0;
deltaH = zeros(1,length(H0));

%% Account for the camera height definition
% This function accomodates camera heights defined with respect to the leak
% location and with respect to the ground. The concentration function is
% written with the camera height defined with respect to the ground.
% The camera height input to the function is stored here, and later
% adjusted to accout for the change of basis if need be.
R03 = cam.R0(3);
detected = zeros(length(Q),1);

N = length(varargin)/2;
if mod(N,1) ~= 0
    display('invalid property definition in LeakDetector');
end
for i=1:2:2*N;
    var = lower(varargin{i});
    if strcmp(var,'heightdefinition')
        if strcmpi(varargin{i+1},'leak')
            deltaH = H0;
        elseif ~strcmp(varargin{i+1}, 'ground')
            display('invalid height definition in LeakDetector');
        end
    elseif strcmp(var,'mindetection')
        minDetection = varargin{i+1};
    else
        display('invalid property name in LeakDetector')
        return
    end
end

%% Determine which leaks are detected
% Loop over each leak
for i = 1:length(Q)
    if Q(i) < minDetection
        detected(i) = 0;
    else
        deltax = (R03 + deltaH(i) - H0(i)) * tan(cam.FOV2/2);
        cam.R0(3) = R03 + deltaH(i);
        cam.R0(2) = y0(i);
        cam.R0(1) = deltax + x0(i);
        
        %% Loop over camera images directly downwind of the leak until either the leak is
        % is detected or the signal begins to decrease.
        iterations = 0;
        while true 
            iterations = iterations + 1;
            if iterations > 5
                display(['Leak Detector Iterations: ' num2str(iterations)])
            end           
            [~, PixelCount] = OverheadPlumeImage(cam, Q(i) , H0(i), x0(i), y0(i),u,Fonethird(i),class);
            if PixelCount == 0 || PixelCountOld > PixelCount || iterations>10
                detected(i) = 0;
                break
            elseif PixelCount > cam.minPixelCount
                detected(i) = 1;
                break
            else
                PixelCountOld = PixelCount;
                cam.R0(1) = cam.R0(1) + deltax;
            end
        end
    end
end