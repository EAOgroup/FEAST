function [ repairedLeaks, remainingLeaks ] = NullLeakFixer( Nleaks, repairRate, deltaT )

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

%Returns the indeces of leaks to be repaired by the null process
% Inputs:
%   Nleaks          the number of leaks
%   repairRate      the number of leaks repaired per leak per day
%   deltaT          length of time of interest (days)
% Ouputs:
%   repairedLeaks   vector of indeces of leaks to be repaired
%   remainingLeaks  vector of indeces of leaks not repaired



% Choose the number of leaks to be repaired by assigining each leak a
% number from a uniform distribution and select leaks with a number above
% the thre threshold set by repairRate*deltaT:
List = rand(1,Nleaks);
repairCount = sum(List<repairRate * deltaT);
% Randomly select the leaks to be repaired
repairedLeaks = randsample(Nleaks,repairCount); %Note: randsample draws without replacement
% Create a vector of remaining leaks.
remainingLeaks = linspace(1,Nleaks,Nleaks);
remainingLeaks(repairedLeaks) = [];
end

