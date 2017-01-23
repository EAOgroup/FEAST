function [leaksFound,leakStruct,repairCost,financeStruct] = DDDetect(DDStruct,financeStruct,processStruct,leakStruct,leaksFound,repairData,timeStruct,k,atmosphereStruct)

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

%% Periodically repair all detected leaks
repairCost = 0;
if mod(timeStruct.time(k),processStruct.repairInterval)<timeStruct.deltaT
    if timeStruct.time(k) ~= 0
        % countFound is a temporary variable to ease notation
        countFound = sum(leakStruct.leaksDetected);
        % account for the cost of locating the leaks
        financeStruct.findCost(k) = financeStruct.locationCost*countFound;
        % calculate the repair cost
        repairCost = sum(datasample(repairData,countFound));
        % update leaksFound and leakStruct
        if countFound>0
            leaksFound = [leaksFound; leakStruct.Q(leakStruct.leaksDetected==1)];
            leakStruct = structfun(@(M)M(leakStruct.leaksDetected==0),leakStruct,'Uniform',0);
        end
    end
    % initialize the list of unfound leaks
    leakStruct.leaksDetected = zeros(length(leakStruct.Q),1);
end
% Define sniffer and leak positioins with respect to the direction of the wind:
x = DDStruct.x*sin(atmosphereStruct.theta) + DDStruct.y*cos(atmosphereStruct.theta);
y = -DDStruct.x*cos(atmosphereStruct.theta) + DDStruct.y*sin(atmosphereStruct.theta);
x0 = leakStruct.x*sin(atmosphereStruct.theta) + leakStruct.y*cos(atmosphereStruct.theta);
y0 = -leakStruct.x*cos(atmosphereStruct.theta) + leakStruct.y*sin(atmosphereStruct.theta);

%% At each time step, iterate through each leak to see if it can be detected
for i = 1:length(leakStruct.Q)
    if leakStruct.leaksDetected(i) == 0
        for m = 1:DDStruct.N
            phi = GaussLeakModel(x(m),y(m),DDStruct.z(m),leakStruct.Q(i),leakStruct.H(i),x0(i),y0(i),atmosphereStruct.u,leakStruct.Fonethird(i),atmosphereStruct.Ry, atmosphereStruct.Rz);
            if phi > DDStruct.Sensitivity
                leakStruct.leaksDetected(i) = 1;
                break
            end
        end
    end
end
