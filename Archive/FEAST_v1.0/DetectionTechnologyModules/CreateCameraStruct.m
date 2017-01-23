function [CameraStruct] = CreateCameraStruct(varargin)

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

%   Creates a struct that defines the relevant properties of a camera
%
% Inputs:   Allowable inputs include:
%       R0 --      1x3 Position of camera [m]
%       n  --      1x3 direction that the camera points
%       top --     1x3 direction of the top of the camera
%       right --   1x3 direction of the right side of the camera
%       npoints -- square root of the number of pixels in the
%                  camera
%       minDetection -- minimum detectable concentration pathlength
%                           [g-m/m^3]
%       FOV1 --     Field of view in the left-right orientation [radians]
%       FOV2 --     Field of view in the top-bottom orientation [radians]
%       theta --    Array of polar angles associated with each pixel
%                   [radians]
%       phi --      Array of azimuth angles associated with each pixel
%                   [radians]
%       minPixelCount -- minimum fraction of pixels with a signal above minDetection
%                        for the leak to be found.
% Outputs:  CameraStruct    	struct with properties set to given inputs
%
% Example usage: CamStruct = createCameraStruct('R0', [5,0,30], 'n', [0,0,-1]);
%
% Author:           Chandler Kemp
% Created:          Apr 8, 2015
% ------------------------------------------------------------------

% Function outline:
%         - Error check the inputs
%         - Create struct base
%         - Save given values into struct
%         - If enough info given (q = 0 or 1), calculate some other properties

% Make sure the right number of inputs are inserted:
if mod(nargin,2) ~= 0  %   if even inputs are required
    error('Even number of inputs required, option-pairs')
end

if nargin > 22 % 2*total num of properties
    error('Too many input arguments in createStreamStruct')
end

%% Create default struct first
R0 = 0;
n = 0;
top = 0;
right = 0;
npoints = 0;
minDetection = 0;
FOV1 = 0;
FOV2 = 0;
theta = 0;
phi = 0;
minPixelCount = 0;

%% Pull variables out of varargin array
for index = 1:2:nargin
    % n will index to the option variables
    switch lower(varargin{index}) % Check what kind of option it is
        case 'r0'
            R0 = varargin{index+1}; % pulls value after 'm'
        case 'n' % lower to be case agnostic
            n = varargin{index+1};
        case 'top' % lower to be case agnostic
            top = varargin{index+1};
        case 'right' % lower to be case agnostic
            right = varargin{index+1};
        case 'npoints'
            npoints = varargin{index+1};
        case 'mindetection'
            minDetection = varargin{index+1};
        case 'fov1'
            FOV1 = varargin{index+1};
        case 'fov2'
            FOV2 = varargin{index+1};
        case 'theta'
            theta = varargin{index+1};
        case 'phi'
            phi = varargin{index+1};
        case 'minpixelcount'
            minPixelCount = varargin{index+1};
        otherwise % Include catch for wrong names
            error(['Invalid Parameter Name: ' (varargin{index})]);
    end
end

%% Make sure camera directions are equal to one and calculate any single missing direction
if max(abs(n))~=0 && max(abs(top))~=0 && max(abs(right)) == 0
    n = n/sqrt(dot(n,n));
    top = top/sqrt(dot(top,top));
    right = cross(n,top);
elseif max(abs(n))==0 && max(abs(top))~=0 && max(abs(right)) ~= 0
    right = right/sqrt(dot(right,right));
    top = top/sqrt(dot(top,top));
    n = cross(right,top);
elseif max(abs(n))~=0 && max(abs(top))==0 && max(abs(right)) ~= 0
    right = right/sqrt(dot(right,right));
    n = n/sqrt(dot(n,n));
    top = cross(right,n);
end

%% Calculate phi and theta if possible and necessary

if npoints ~= 0 && FOV1 ~= 0 && FOV2 ~= 0 && theta == 0;
    [theta, phi] = PixelMaker(npoints,FOV1,FOV2);
elseif theta ~= 0
    if size(theta) ~= [npoints,npoints]
        error('size(theta) is not [npoints,npoints]');
    end
elseif phi ~= 0
    if size(phi) ~= [npoints,npoints]
        error('size(phi) is not [npoints,npoints]');
    end
end


%% Save away final values
CameraStruct.R0 = R0;
CameraStruct.n = n;
CameraStruct.top = top;
CameraStruct.right = right;
CameraStruct.npoints = npoints;
CameraStruct.minDetection = minDetection;
CameraStruct.FOV1 = FOV1;
CameraStruct.FOV2 = FOV2;
CameraStruct.theta = theta;
CameraStruct.phi = phi;
CameraStruct.minPixelCount = minPixelCount;
CameraStruct.imageWidth = tan(FOV1/2)*R0(3)*2;

end % end function