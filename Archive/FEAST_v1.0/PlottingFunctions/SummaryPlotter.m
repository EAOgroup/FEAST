function [ varargout ] = SummaryPlotter( Directory )

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

% Plots summary of results in the folder 'Directory'
% Inputs
%   Directory     Path to a folder to be analyzed. 'Directory' must include
%                 results from atleast one simulation
% Outputs: No outputs are required. Optional outputs in order are:
%   Total NPV           Struct of total marginal NPV [k$/well]
%   Standard Deviation  Struct of standard deviation in NPV [k$/well]
%   Repair              Struct of repair costs [k$/well]

ColorSet = [ 140,21,21; 0,152,219 ; 0,155,118;  178,111,22; 83,40,79;  0,0,0; 77,79,83];
ColorSet = ColorSet/max(max(ColorSet));
set(0,'DefaultAxesColorOrder',ColorSet)

% Load example file for settings
files = dir(Directory);
load([Directory '/' files(3).name])

%Calculate NPV info
[~, nullNPV,~, timeSeries] = ResultsAnalysis(Directory);
Nruns = length(nullNPV);
Ntech = length(timeSeries(:,1,1))-1;

%% Marginal NPV calcs

NPV = nullNPV;
Capital = zeros(Ntech,1); Maintenance = zeros(Ntech,1); Repair = zeros(Ntech,1); Finding = zeros(Ntech,1); Leakage = zeros(Ntech,1); Total = zeros(Ntech,Nruns);
% Calculate the average NPV of each component
for i = 1:Nruns
    Capital = Capital + NPV(i).Capital/Nruns;
    Maintenance = Maintenance +NPV(i).Maintenance/Nruns;
    Repair = Repair + NPV(i).Repair/Nruns;
    Finding = Finding + NPV(i).Finding/Nruns;
    Leakage = Leakage + NPV(i).savedLeakage/Nruns;
    Total(:,i) = NPV(i).Total;
end

standardDev = std(Total,0,2);
Total = mean(Total,2);

for i = 1:length(legendString)
    legendString2{2*i-1} = {' '};
    legendString2{2*i} = legendString{i};
end

%% Marginal stacked bar
figure(1)
clf
% Plot the costs
H = bar(1:Ntechnologies-1,[Capital(1:Ntechnologies-1), Finding(1:Ntechnologies-1), Repair(1:Ntechnologies-1), Maintenance(1:Ntechnologies-1)],'stack');
set(gca,'XTickLabel',legendString(1:Ntechnologies-1))
hold all
% Add the benefit bars
H(end+1) = bar(Leakage(1:Ntechnologies-1),'stack','c');
% Customize the color of the bars
P=findobj(gca,'type','patch');
for n=1:length(H)
    set(H(n),'facecolor',ColorSet(n,:))
end

%Add axis labels, legend, etc.
set(gca,'XLim',[0.5 Ntechnologies-1+0.5])
H(end+1) = errorbar(Total(1:Ntechnologies-1),standardDev(1:Ntechnologies-1)/sqrt(Nruns),'.','MarkerSize',12,'Color',ColorSet(6,:));
ylabel('Cost and benefit present value ($1k/well)')
xlabel('Method')
plotfixer

legend1 = legend(H,'Capital Cost','Finding Cost','Repair Cost','Maintenance Cost','Value of gas saved','NPV');
set(legend1,...
    'Position',[0.594137286302137 0.152957280289578 0.353067814854683 0.186304128902316],...
    'LineWidth',1.5);
