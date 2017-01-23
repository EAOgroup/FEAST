function [ noRepairNPV, nullNPV, leaksFoundCum, timeSeries ] = ResultsAnalysis(path)

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

% This function reads results saved in .mat files from a
% directory and summarizes results.
%
% Inputs:
% directory     path to a directory full of results run with the same
%               Ntechnologies, timeSteps, and SiteCount
% Outputs;
% noRepairNPV   Array of structs of no repair NPVs for every technology in every
%                   result file [$1k/well]
% nullNPV       Array of structs of null NPVs for every technology in every
%                   result file [$1k/well]
% leaksFound    List of leaks found using each technology [g/s]
% timeSeries    Array with dimension (Ntechnologies + 1, timeSteps, Nresuls) storing
%               the leakage time series for every technology in every
%               result file [g/s]

files = dir(path);
load([path '/' files(3).name ])
timeSeries = zeros(Ntechnologies+1,timeStruct.timeSteps,length(files)-2);

for j = 1:Ntechnologies
    if strcmpi(legendString{j},'Null')
        leaksFoundCum.(legendString{j}) = nullRepaired(j).Q;
    else
        leaksFoundCum.(legendString{j}) = leaksFound(j).Q;
    end
    [noRepairNPV(1),nullNPV(1)] = NPVCalculator([path '/' files(3).name]);
    timeSeries(:,:,1) = Leakage;
end

for iindex = 4:length(files)
    load([path '/' files(iindex).name]);
    timeSeries(:,:,iindex-2) = Leakage;
    for j = 1:Ntechnologies
        leaksFoundCum.(legendString{j}) = [leaksFoundCum.(legendString{j}); leaksFound(j).Q];
    end
    % Calculate capital, operating, and repair NPV
    % First initialize variables
    [noRepairNPV(iindex-2), nullNPV(iindex-2)] = NPVCalculator([path '/' files(iindex).name]);
end