function [leaksFound,leakStruct,repairCost,financeStruct] = MIRDetect(MIRStruct,financeStruct,processStruct,leakStruct,leaksFound,repairData,timeStruct,k,atmosphereStruct)

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

% MIRDetect determines which leaks can be identified in the MIR program.
% Inputs:

repairCost = 0;
if mod(timeStruct.time(k),processStruct.surveyInterval)<timeStruct.deltaT
    % Sort leaks by leak size
    % (Note that in the MIR module, the leak location does not affect
    % detection. Therefore leak detection under a given set of atmospheric
    % conditions depends on the leak flux exclusively.
    [~, SortVec] = sort(leakStruct.Q,'descend');
    leakStruct = structfun(@(M)M(SortVec),leakStruct,'Uniform',0);
    detected = zeros(length(leakStruct.Q),1);
        % determine which leaks are detected:
    for i = 1:length(leakStruct.Q)
        leak = structfun(@(M)M(i),leakStruct,'Uniform',0);
        detected(i) = LeakDetector(MIRStruct,leak,atmosphereStruct,'HeightDefinition','leak','minDetection',MIRStruct.Qmin);
        if detected(i) == 0
            break
        end
    end
    countFound = sum(detected);
    
    if countFound > 0
        repairCost = sum(datasample(repairData,countFound));
        leaksFound = [leaksFound; leakStruct.Q(detected==1)];
        % Only save leaks that were not found
        leakStruct = structfun(@(M)M(detected==0),leakStruct,'Uniform',0);
    end
end