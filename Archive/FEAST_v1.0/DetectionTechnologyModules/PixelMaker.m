function [ theta phi]=PixelMaker(npoints,FOV1,FOV2)

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

%Defining the grid of pixels in an IR camera. I define a cartesian grid here
% to get a useful  set of points and later convert to spherical coordinates.
% The units of the grid end up being radians. The result is something akin 
% to pixels: each grid location is represented as a square in 2D, but it 
% represents an integral along a path at a specific angle

% Inputs:
%   npoints     square root of the total number of pixels
%   FOV1        horizontal (relative to the camera's "top" vector) field of
%               view angle [radians]
%   FOV2        vertical (relative to the camera's "top" vector) field of
%               view angle [radians]         

% Outputs:
%   theta       array of pixel polar angles in spherical coordinates
%               [radians]
%   phi         array of pixel azimuth angles in spherical coordinates
%               [radians]


x=linspace(-FOV1/2,FOV1/2,npoints);
y=linspace(-FOV2/2,FOV2/2,npoints);
phi=zeros(npoints);
theta=zeros(npoints);
dx=x(2)-x(1);
dy=y(2)-y(1);

for i = 1:npoints %iterate over the grid in x
    for j = 1:npoints %iterate over the grid in y
        %theta defines the polar angle from the camera to the current point
        theta(j,i)=sqrt(x(i)^2+y(npoints+1-j)^2);
        %This if statement takes care of the multivalued property of
        %inverse tangent which is not captured by atan
        if x(i)<0
            phi(j,i)=atan(y(npoints+1-j)/x(i))+ pi;
        else
            phi(j,i)=atan(y(npoints+1-j)/x(i));
        end
    end
end