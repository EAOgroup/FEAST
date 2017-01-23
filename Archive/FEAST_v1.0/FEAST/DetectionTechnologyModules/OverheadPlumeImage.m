function [ Image, pixelCount  ] = OverheadPlumeImage( cam, Q , H0, x0, y0,u,Fonethird,class )

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

%OverheadPlumeImage calculates the concentration pathlength associated with
%each pixel in a simulated camera and uses it to determine the fraction of
%pixels that can detect the leak. The camera must be looking straight down
%at a leak from above.

% Inputs:
%   Cam -- Struct of camera properties. Required fields include:
        % R0     3D position of the camera [m]
        % n      direction that the camera points
        % top    direction that the top of the camera points
        % CamRight direction that the right side of the camera points
        % npoints square root of the number of pixels to be simulated
        % minDetection minimum detectable concentration pathlength [g-m/m^3]
        % theta  array of polar angles associated with each pixel in the camera
        %        [radians]
        % phi    array of azimuth angles associated with each pixel in the
        %        camera [radians]
%   Non-struct inputs:
    % Q      vector containing the flux of each leak [g/s]
    % H0     height of the leak [m], 
    % x0, y0 define the location of each leak [m]
    % u      windspeed [m/s]
    % Fonethird Buoyancy flux raised to the one-third power [(m^4/s^3)^1/3]
    % class  Defines the stability class of the atmosphere surrounding the
    %        leak.
    
% Outputs:
    % Image    array storing the concentration pathlength of methane for
    %               each pixel [g-m/m^3]
    % PixelCount fraction of pixels with a concentration pathlength greater
    %               than the minimum detecatable value.
    
%% Proceed with the calculation
% Define integration steps and pull common variables from the cam struct
nsteps = 15; 
npoints = cam.npoints; R0 = cam.R0;

% Initialize variables
Image=zeros(npoints,npoints);
z0=H0;
%sigmay and sigmaz are estimates of the standard deviation of the plume at
%the location along the path of each pixel where the concentration is at
%its maximum.
sigmay = zeros(npoints,npoints);
sigmaz = zeros(npoints,npoints);

%The following five parameters will be used to calculate sigmaz and
%sigmay as a function of x within GaussLeakModel. From Seinfeld page 866.

Ry = [0.23 0.18 0.11 0.08 0.058 0.04];
Rz = [0.18 0.12 0.07 0.05 0.035 0.025];
Ry = Ry(class); Rz = Rz(class);

%s will store an estimate of the distance from a pixel to the point
%along the path imaged by that pixel with the highest concentration.
s = zeros(npoints,npoints);

% Define the rate of change of x, y, z along the path associated with each
% pixel.
dx=sin(cam.theta).*(cam.right(1)*cos(cam.phi));
dy=sin(cam.theta).*(cam.top(2)*sin(cam.phi));
dz=cam.n(3)*cos(cam.theta);


%z0 gives the centerline height of the plume at the location of the camera.
if R0(1)-x0>=0
    z0 = z0+1.6 * Fonethird * (R0(1)+(R0(3)-H0)*dx-x0).^(2/3)./u;
end
%R01 stores estimates of the x position of the path of each pixel at the
%point of maximum concentration along the path.
R01 = R0(1)+dx.*(R0(3)-z0);
%pL stores the indices of the points with R01>0
pL = find(R01-x0>0);

sigmay(pL) = Ry*(R01(pL)-x0);
sigmaz(pL) = Rz*(R01(pL)-x0);

% Choose sigmaz for stepsize since sigmaz<sigmay
stepsize = 8*sigmaz/nsteps;
%This gives an estimate of the point of max concentration.
s(pL) = -(sigmaz(pL).^2*R0(2).*dy(pL)+dz(pL).*sigmay(pL).^2.*(R0(3)-z0(pL)))./(sigmaz(pL).^2.*dy(pL).^2+sigmay(pL).^2.*dz(pL).^2);
s(dy.^2+dz.^2 < .0001)=R0(3)-H0;
s = max(s-nsteps/2*stepsize,0);
x = R0(1) + dx.*s; y=R0(2)+dy.*s; z = R0(3) + dz.*s;
% Calculate the concentration at points near the maximum concentration to
% estimated the total concentration pathlength.
for i=1:nsteps
    x = x+dx.*stepsize; y = y+dy.*stepsize; z=z+dz.*stepsize;
    Conc = GaussLeakModel(x,y,z,Q,H0,x0,y0,u,Fonethird,Ry, Rz);
    Image = Image + Conc;
    xvec(i) = x(npoints/2,npoints/2+1);
    yvec(i) = y(npoints/2,npoints/2+1);
    zvec(i) = z(npoints/2,npoints/2+1);
end
Image=Image.*stepsize;

pixelCount = sum(sum(Image>cam.minDetection))/npoints^2;