function [ Phi ] = GaussLeakModel(x,y,z,Q,H0,x0,y0,u,Fonethird,Ry, Rz)

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

%GaussLeakModel returns the concentration of a pollutant dispersing from one
%or more point sources. 

% Inputs:
    % x,y,z  array of points at which to calculate the concentration [m]
    % Q      vector containing the flux of each leak [g/s]
    % H0     vector defining the height of the leaks [m], 
    % x0, y0 Location of each leak [m]
    % u      windspeed [m/s]
    % Fonethird     Buoyancy flux raised to the one-third power [(m^4/s^3)^1/3]
    % Ry, Rz        scalars defining the rate of change of sigmay and sigmaz
    
% Outputs:
    % Phi    array storing the concentration of methane at the desired
    %        coordinates [g/m^3]

%% Definitions
%sigmay is the standard deviation of the concentration in the y direction.
%sigmaz is the standard deviation of the concentration in the z direction.

Phi=zeros(size(x)); % Phi is the concentration at each location specified by x, y and z.

%%Model
%sigmay and sigmaz are calculated as a function of x based on the Pasquill
%stability category. They are linear fits to the 100 meter data on the
%Pasquill Gifford curves.

% Iterate through each leak
for i=1:length(Q)
    normx=x-x0(i);
    % poslocations is used to only do calculations in the region that is
    % not known to be zero for all leaks
    poslocations = find(z > 0 & normx > 0);
    if ~isempty(poslocations)
        normxpos=normx(poslocations);
        %sigmay and sigmaz are arrays of dimension normxpos
        sigmay = Ry*normxpos;
        sigmaz = Rz*normxpos;
        % a and b are variables used in later calculations to account for
        % the buoancy of the plume.
        b=-H0(i)-1.6 * Fonethird * normxpos.^(2/3)/u;
        a=z(poslocations)+b;
        %Concentration Calculation
        fy = exp(-(y(poslocations)-y0(i)).^2./(2*sigmay.^2));
        fz = exp(-a.^2./(2*sigmaz.^2));
        fg = exp(-(a-2*b).^2./(2*sigmaz.^2));
        Phi(poslocations) = Phi(poslocations)+Q(i)./(u*2*pi.*sigmaz.*sigmay).*fy.*(fz+fg); %g/m^3
    end
end
end

