function [ noRepairNPV, nullNPV ] = NPVCalculator( filePath )

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

%This function opens a results file and calculates the NPV components
%   Inputs:
%           filePath    Path to a results file
%   Outputs (all units 1k$/site):
%       noRepairNPV is a struct with the following fields [$1k/well]:
%           Capital     NPV of capital investments
%           Finding     NPV of operating costs
%           Maintenance NPV of maintenance costs
%           Repair      NPV of repair costs
%           savedLeakage NPV of savedLeakage
%           total       Total NPV
%       nullNPV has the same structure as NPV but stores the null NPV
%           rather than the noRepairNPV

load(filePath)
SiteCount = GasField.SiteCount;
NullIndex = find(strcmpi('Null',legendString));

discountArray = ones(Ntechnologies,1)*(1+econStruct.discountRate).^(timeStruct.time/365);

for i = 1:Ntechnologies
    savedLeakageValuePV(i) = sum((Leakage(end,:)-Leakage(i,:))./discountArray(1,:))*timeStruct.deltaT*24*3600*econStruct.gasPrice;
    savedLeakageValuePVnull(i) = sum((Leakage(NullIndex,:)-Leakage(i,:))./discountArray(1,:))*timeStruct.deltaT*24*3600*econStruct.gasPrice;
end

findCostPV = Finding./discountArray;
repairCostPV = repairCost./discountArray;
maintenanceCostPV = Maintenance./discountArray;
capitalCostPV = Capital./discountArray;
repairCostPVnull = (repairCost-ones(Ntechnologies,1)*repairCost(NullIndex,:))./discountArray;

noRepairNPV.Capital = -sum(capitalCostPV,2)/SiteCount/1000;
noRepairNPV.Finding = -sum(findCostPV,2)/SiteCount/1000;
noRepairNPV.Repair = -sum(repairCostPV,2)/SiteCount/1000;
noRepairNPV.Maintenance = -sum(maintenanceCostPV,2)/SiteCount/1000;
noRepairNPV.savedLeakage = savedLeakageValuePV'/SiteCount/1000;
noRepairNPV.Total  = noRepairNPV.Capital + noRepairNPV.Finding + noRepairNPV.Repair + noRepairNPV.Maintenance + noRepairNPV.savedLeakage;
nullNPV.Capital = noRepairNPV.Capital;
nullNPV.Finding = noRepairNPV.Finding;
nullNPV.Repair = -sum(repairCostPVnull,2)/SiteCount/1000;
nullNPV.Maintenance = noRepairNPV.Maintenance;
nullNPV.savedLeakage = savedLeakageValuePVnull'/SiteCount/1000;
nullNPV.Total  = nullNPV.Capital + nullNPV.Finding + nullNPV.Repair + nullNPV.Maintenance + nullNPV.savedLeakage;

end


