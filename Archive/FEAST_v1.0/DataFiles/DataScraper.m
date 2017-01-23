%% This script includes the code to produce all of the 'scraped' datasets.
%% Fort Worth Scraper:

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

%this is the script used to produce 'FortWorthScraped' from 'fortworth'
%In 'FortWorthScraped', all of the methane fluxes reported as zero are
%scraped from the dataset. All flux data from records with a reported hiflow
%methane concentration of .001 are also scraped from the dataset (.001 was used  
%as a default value for leaks with no measurable methane concentration in the hiflow
%sampler exhaust).
%Note that The IR camera was used to survey all
%components, while the traditional Method 21 survey was only conducted on
%10% of all components.

% First delete rows with percentCFM == 0.001
data=importdata('FortWorthPercentCFMCameraBoolPPM.mat');
data(data(:,1) == 0.001,:) = [];
Nwells = 1138;

% Now delete rows with percentCFM == 0;
data(data(:,1) == 0,:) = [];

% convert units from CFM to g/s (assuming standard conditions)
data(:,1) = data(:,1) * 0.0283/60*1e5/8.314/293*16;

% Now make a vector of leaks found with and without the IR camera
FortWorthScraped.IR = data(data(:,2) == 0,1);
FortWorthScraped.FID = data(data(:,2) == 1,1);
FortWorthScraped.allLeaks = data(:,1);

% Save the vector of measured concentrations
FortWorthScraped.concentrations = data(:,3);

% Also save the rate at which leaks were found by the IR camera
% [leaks/well]
FortWorthScraped.camLeakRate = length(FortWorthScraped.IR)/Nwells;

% and the rate at which leaks werefound generally (accounting for the fact
% that while the camera was used to find leaks throughout every facility,
% the other leaks were found while searching just ten percent of the
% facilities.

FortWorthScraped.genLeakRate = (length(FortWorthScraped.FID)+length(FortWorthScraped.IR)/10)/(Nwells/10);
FortWorthScraped.avgLeak = 0.9*mean(FortWorthScraped.FID) + 0.1*mean(FortWorthScraped.IR);
    

%% Allen Scraper:
%Leaks with a reported flux of zero are scraped from the dataset.
% AllenSCFMMethaneScraped=importdata('AllenSCFMMethane.mat');
% i=1;
% while i <= length(AllenSCFMMethaneScraped)
% if AllenSCFMMethaneScraped(i)<=0
%     AllenSCFMMethaneScraped(i)=[];
% else
%     i=i+1;
% end
% end
