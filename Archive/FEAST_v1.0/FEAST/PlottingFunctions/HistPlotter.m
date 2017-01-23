function [  ] = HistPlotter( directory )

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

% Plots histogram of leaks found by each technology based on results in the folder 'Directory'
% Inputs
    % Directory     Path to a folder to be analyzed
    
% Outputs
%    None
%% Define color scheme
ColorSet = [ 140,21,21; 0,152,219 ; 0,155,118;  178,111,22; 83,40,79;  0,0,0; 7\7,79,83];
ColorSet = ColorSet/max(max(ColorSet));
set(0,'DefaultAxesColorOrder',ColorSet)

[~,~,Found,~] = ResultsAnalysis(directory);
% Load sample file for GasField settings
files = dir(directory);
load([directory '/' files(3).name])

% Create a fit
nbins = 30;
ParamsLeakDist = fitdist(GasField.LeakData,'Lognormal');
[counts, centers]=hist(GasField.LeakData,nbins);
siteCount = (length(files)-2)*GasField.SiteCount;
delta = centers(2)-centers(1);
x=linspace(min(GasField.LeakData),max(GasField.LeakData)-delta/2,2*nbins);
y = lognpdf(x,ParamsLeakDist.mu,ParamsLeakDist.sigma)*delta*(length(GasField.LeakData)+timeStruct.endTime*GasField.LeakProductionRate*GasField.ComponentCount)*(length(files)-2);

%Create a histogram for the first four LDAR technologies. May adjust field names for
%custom technology modules.

names = fields(Found);
for i =1:length(names)
figure(i)
clf
[count, centers] = hist(Found.(names{i})(Found.(names{i})~=0),nbins);
hb1 =  bar(centers,count/siteCount,1,'FaceColor',ColorSet(i,:));
% title('Leaks found with AIR')
xlabel('Leak size [g/s]')
ylabel('Leaks found per well')
hold on
hb2 = plot(x,y/siteCount,'Color',ColorSet(6,:));
legend([hb1,hb2],'Histogram of Leaks Found','Fit to Leak Distribution')
ylim([0,2])
text(0.3,1.15,names{i})
plotfixer
end

