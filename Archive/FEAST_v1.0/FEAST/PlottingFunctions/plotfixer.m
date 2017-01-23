%Plot Fixer
%Originally Written by: Matt Svrcek  12/05/2001
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


%Run this script after generating the raw plots.  It will find
%all open figures and adjust line sizes and text properties.

%Change the following values to suit your preferences.  The variable
%names and comments that follow explain what each does and their options.

plotlsize = 2; %thickness of plotted lines, in points
axislsize = 1.5; %thickness of tick marks and borders, in points
markersize = 8;  %size of line markers, default is 6

%font names below must exactly match your system's font names
%check the list in the figure pull down menu under Tools->Text Properties
%note, the script editor does not have all the fonts, so use the figure menu

axisfont = 'Helvetica'; %changes appearance of axis numbers
axisfontsize = 18;            %in points
axisfontweight = 'normal';    %options are 'light' 'normal' 'demi' 'bold' 
axisfontitalics = 'normal';   %options are 'normal' 'italic' 'oblique'

legendfont = 'Helvetica'; %changes text in the legend
legendfontsize = 18;
legendfontweight = 'normal';
legendfontitalics = 'normal';

labelfont = 'Helvetica';  %changes x, y, and z axis labels
labelfontsize = 18;  
labelfontweight = 'normal'; 
labelfontitalics = 'normal';

titlefont = 'Helvetica';  %changes title
titlefontsize = 18;
titlefontweight = 'normal';
titlefontitalics = 'normal';

textfont = 'Helvetica';   %changes text
textfontsize = 18;
textfontweight = 'normal';
textfontitalics = 'normal';


%stop changing things below this line
%----------------------------------------------------
axesh = findobj('Type', 'axes');
legendh = findobj('Tag', 'legend');
lineh = findobj(axesh, 'Type', 'line');
axestexth = findobj(axesh, 'Type', 'text');

set(lineh, 'LineWidth', plotlsize)
set(lineh, 'MarkerSize', markersize)
set(axesh, 'LineWidth', axislsize)
set(axesh, 'FontName', axisfont)
set(axesh, 'FontSize', axisfontsize)
set(axesh, 'FontWeight', axisfontweight)
set(axesh, 'FontAngle', axisfontitalics)
set(axestexth, 'FontName', textfont)
set(axestexth, 'FontSize', textfontsize)
set(axestexth, 'FontWeight', textfontweight)
set(axesh, 'Box','on')
for(i = 1:1:size(axesh))
   legend(axesh(i))
   set(get(gca,'XLabel'), 'FontName', labelfont)
   set(get(gca,'XLabel'), 'FontSize', labelfontsize)
   set(get(gca,'XLabel'), 'FontWeight', labelfontweight)
   set(get(gca,'XLabel'), 'FontAngle', labelfontitalics)
   set(get(gca,'YLabel'), 'FontName', labelfont)
   set(get(gca,'YLabel'), 'FontSize', labelfontsize)
   set(get(gca,'YLabel'), 'FontWeight', labelfontweight)
   set(get(gca,'YLabel'), 'FontAngle', labelfontitalics)
   set(get(gca,'ZLabel'), 'FontName', labelfont)
   set(get(gca,'ZLabel'), 'FontSize', labelfontsize)
   set(get(gca,'ZLabel'), 'FontWeight', labelfontweight)
   set(get(gca,'ZLabel'), 'FontAngle', labelfontitalics)
   set(get(gca,'Title'), 'FontName', titlefont)
   set(get(gca,'Title'), 'FontSize', titlefontsize)
   set(get(gca,'Title'), 'FontWeight', titlefontweight)
   set(get(gca,'Title'), 'FontAngle', titlefontitalics)
end
