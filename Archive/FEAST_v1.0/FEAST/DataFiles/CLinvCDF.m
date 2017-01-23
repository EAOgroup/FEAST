function [ x ] = CLinvCDF( y )

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

%This function takes a vector of numbers between 0 and 1 and returns an
%equal length vector of values calculated from the input using an inverse
%cumulative distribution function (CDF).
% The inverse CDF function was calculated based on the mean, median, min, and max
% leak repair costs provided by Carbon Limits for leak repair.
% Inputs:
%   y   Vector of values between 0 and 1 representing the position in a
%           cumulative distribution
% Outputs:
%   x   Vector of values computed by function

if min(y) < 0 || max(y)>1
    error('Input values must be between zero and 1')
end

% Carbon Limits data for leak repair costs:
minVec = [20 15 20 20];
medVec = [50 50 125 50];
maxVec = [5500 5000 1000 2000];
meanVec = [90 56 189 129];
componentRatio = [10575 23577 1081 1106];
componentRatio = componentRatio/sum(componentRatio);

% Calculate the weighted average of each piece of data
minim = sum(minVec.*componentRatio);
median = sum(medVec.*componentRatio);
maxim = sum(maxVec.*componentRatio);
meanVal = sum(meanVec.*componentRatio);
a = minim; b = median; c = a + 2*(b-a); d = maxim;

% Create the inverse CDF. 
% The CDF is computed by assuming a probability density function (PDF) that is a
% modified triangular distribution with a long tail for rare, high cost
% repairs. The PDF roughly has this shape: /\_
syms h1 h2
eqn1 = meanVal == h1/(b-a)*(b^3/3-a*b^2) + h1/2*(c^2-b^2) - (h1-h2)/(c-b)*(c^3/3-b*c^2/2) + h2/2*(d^2-c^2) - h2/(d-c) * (d^3/3-c*d^2/2);
eqn2 = 1 == (b-a)*h1/2 + (c-b)*(h1-h2)/2 + (c-b)*h2 +(d-c)*h2/2;
S = solve([eqn1,eqn2],h1,h2);

% h1 and h2 are the height of the peak and the beginning of the tail in the
% PDF.
h1 = double(S.h1);
h2 = double(S.h2);
x = nan(size(y));

%C1 and C2 correspond to the peak of the PDF and the beginning of the tail,
%respectively.
C1 = h1/(b-a)*((b^2/2-a*b)+a^2/2);
C2 = C1 + (h1*c-(h1-h2)/(c-b)*(c^2/2 - b*c)) - h1*b-(h1-h2)/(c-b)*b^2/2;

C1 = double(C1); C2 = double(C2);
theta = h1/(b-a);
phi = (h1-h2)/(c-b);
xi = h2/(d-c);

% The x values are computed based on the slope of each section of the pdf.
x(y<C1) = (a*theta + sqrt(a^2*theta^2 -4*theta/2*(theta/2*a^2-y(y<C1))))/theta;
x(y>=C1 & y < C2) = (-(h1+b*phi) + sqrt((h1+b*phi)^2 + 4*phi/2*(-h1*b-phi*b^2/2 - y(y>=C1 & y < C2)+C1)))/(-phi);
x(y >= C2) = (-(h2+c*xi) + sqrt((h2+c*xi)^2 + 4*xi/2*(-h2*c-xi*c^2/2 - y(y>=C2)+C2)))/(-xi);

