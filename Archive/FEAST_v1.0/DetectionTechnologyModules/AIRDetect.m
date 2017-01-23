function [leaksFound,leakStruct,repairCost,financeStruct] = AIRDetect(AIRStruct,financeStruct,processStruct,leakStruct,leaksFound,repairData,timeStruct,k,atmosphereStruct)

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

repairCost = 0;
if mod(timeStruct.time(k),processStruct.surveyInterval)<timeStruct.deltaT
    % determine which leaks are detected:
    detected = LeakDetector(AIRStruct,leakStruct,atmosphereStruct,'minDetection',AIRStruct.Qmin);
    countFound = sum(detected);    
    if countFound > 0
        repairCost = sum(datasample(repairData,countFound));
        leaksFound = [leaksFound; leakStruct.Q(detected==1)];
        % Only save leaks that were not found
        leakStruct = structfun(@(M)M(detected==0),leakStruct,'Uniform',0);
    else
        repairCost = 0;
    end
end