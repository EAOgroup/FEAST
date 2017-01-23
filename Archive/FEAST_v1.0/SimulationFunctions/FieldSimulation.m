function [] = FieldSimulation(varargin)

%-------------------------------------------------------------------------------------------------------------------
% Fugitive Emissions Abatement Simulation Toolkit aka FEAST
% Copyright {2016} {Chandler E. Kemp; Arvind P. Ravikumar; Adam R. Brandt} 

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.
%--------------------------------------------------------------------------------------------------------------------

% This script runs a natural gas field simulation. It calculates the
% leakage at a series of timesteps under a number of LDAR programs
% Inputs:
%   varargin may provide paths to two optional input files:
%       1.) A path to a series of generated leaks
%       2.) A path to an initial realization of a Gas Field
%% Outputs: These variables are saved to a file at the end of the function
%   Leakage         time series of leakage [g/s]
%   econStruct     struct of economic variables with the following fields:
%       gasPrice        cost of gas [$/g]
%       discountRate    real discount rate for NPV analysis [fraction per
%       year]
%   atmStruct       Atmospheric data at every timestep including:
%       u               windspeed [m/s]
%       theta           wind direction [radians from North]
%       class           Atmospheric stability class
%       Ry              Constant determining sigma_y in the Gaussian plume
%                       model (depends on stability class)
%       Rz              Constant determining sigma_z in the Gaussian plume
%                       model (depends on stability class)
%   leakStruct      Array of structs defining leaks existing at the end of
%                   the simulation for each technology module. Has the
%                   following fields:
%       Q               Size of leaks [g/s]
%       x               x position of leaks [m]
%       y               y position of leaks [m]
%       H               Vertical position of leaks [m]
%       Fonethird       Buoyancy flux parameter raised to the one-third
%                           for every leak [m^{4/3}/s]
%       leaksDetected   Leaks detected but not yet repaired [binary]
%   newLeaks        Array of structs defining leaks created at each time 
%                   step. Has the following fields:
%       Q               Size of leaks [g/s]
%       x               x position of leaks [m]
%       y               y position of leaks [m]
%       H               Vertical position of leaks [m]
%       Fonethird       Buoyancy flux parameter raised to the one-third
%                           for every leak [m^{4/3}/s]
%       leaksDetected   Leaks detected but not yet repaired [binary]
%   leaksFound      Array of structs leaks that were found and repaired by
%                   each LDAR program. The structs have one field:
%       Q               Flux from each repaired leak [g/s]
%   nullRepaired    Array of structs of leaks that were found and repaired
%                   by the null process in each LDAR program. Each struct
%                   has one field: 
%       Q               Flux from each repaired leak [g/s]
%   processStruct    Array of structs defining process variables for each
%                   detection technology. Fields depend on the type of technology.
%   techStruct      Array of structs defining properties of each detection
%                   technology. Fields depend on the technology.
%   timeStruct      Struct of time variables for the simulation
%       endTime         Time at the end of the simulation [days]
%       timeSteps       Number of time steps in the simulation
%       time            Array of times considered in the simulation [days]
%       deltaT          Length of one time step [days]
%   Capital         Array of capital expenditures. First index determines 
%                   the technology, second index determines the timestep.
%   Finding         Array of finding costs. First index determines 
%                   the technology, second index determines the timestep.
%   Leakage         Array of leakage. First index determines 
%                   the technology, second index determines the timestep.
%                   Includes leakage in the "No Repair" scenario in the
%                   final row.
%   Maintenance     Array of maintenance costs. First index determines 
%                   the technology, second index determines the timestep.
%   Ntechnologies   Number of distinct LDAR programs in the simulation.
%   fileName        Path where results were originally saved.
%   leakList        List of leaks generated in the simulation (including 
%                   initial leaks) [g/s]
%   legendString    Cell array of LDAR program names
%   repairCost      Array of repair costs. First index determines the
%                   technology, second index determines the timestep.

%% Set up
% Defines a new seed for random variables.
rng('shuffle')

% Add necessary folders
addpath([pwd '/DataFiles'])
addpath([pwd '/DetectionTechnologyModules'])
addpath([pwd '/InitializationFunctions'])
addpath([pwd '/SimulationFunctions'])

%% Simulation parameters
endTime = 10*365; % days
timeSteps = 4000;
deltaT = endTime/(timeSteps - 1); %days
time = linspace(0,endTime,timeSteps);
timeStruct.endTime = endTime; timeStruct.timeSteps = timeSteps; timeStruct.time = time; timeStruct.deltaT = deltaT;

% Initiate variables
Ntechnologies = 5;
Leakage = zeros(Ntechnologies+1,timeSteps);
legendString = {'AIR','DD','MIR','FID','Null'};

%% Initialize variables
for i = 1:2:length(varargin)-1
    argin = lower(varargin{i});
    if strcmp(argin,'gasfield')
        GasField = varargin{i+1};
    elseif strcmp(argin,'leakseries')
        leakSeries = varargin{i+1};
    elseif strcmp(argin,'atmosphere')
        atmStruct = varargin{i+1};
    end
end

if ~exist('GasField','var')
    % Create a gas field
    GasField = GasFieldInitiator();
end
% Initiate the list of leaks
leakList = GasField.LeakSample;

% Load the economic parameters (econStruct is not used in the simulation,
% but is important to interpret results once the simulation is complete)
econStruct = FinancialAssumptions();

% Load repair cost data (Either Fernandez or Carbon Limits data)
repairdata = importdata([pwd '/DataFiles/FernandezRepairCost.mat']);
%repairdata = CLinvCDF(rand(1e4,1));

% Create atmosphere data
if ~exist('atmStruct','var')
    atmStruct = AtmosphereInitiator(timeSteps);
end

% Initiate the vector of repair costs
repairCost = zeros(Ntechnologies,timeSteps);

%% Technology initialization
% Loop through each technology
for i = 1:Ntechnologies
    % Define a function handle for the current technology:
    fh = str2func(legendString{i});
    % Call the function to intialize the current technology
    [techStruct{i}, financeStruct{i}, processStruct{i}] = fh(GasField, timeStruct);
    leakStruct(i).Q = GasField.LeakSample;        % List of current leaks [g/s]
    % List of leak positions [m]
    leakStruct(i).x = GasField.x0; leakStruct(i).y = GasField.y0; leakStruct(i).H = GasField.H0;
    % Leak buoyancy flux raised to the one-third power
    leakStruct(i).Fonethird = GasField.Fonethird;
    % LeaksFound [g/s]
    leaksFound(i).Q = [];
    % List of leaks detected but not yet repaired
    leakStruct(i).leaksDetected = zeros(length(leakStruct(i).Q),1);
    % LeaksFound with the Null process
    nullRepaired(i).Q = [];
end

%% Simulation
for k = 1:length(time)
    if mod(k,100) == 0
        display(['Currently evaluating timestep number ' num2str(k)])
    end
    % Generate new leaks
    if exist('leakSeries','var') % Load leaks from input data if possible
        newLeaks(k) = leakSeries(k);
    else
        newLeaks(k) = LeakCreater(GasField.ComponentCount, deltaT, GasField.LeakProductionRate, GasField.LeakData);
    end
    % Set variables for the timestep
    atmosphere = structfun(@(M)M(k),atmStruct,'Uniform',0);
    %% No repair scenario
    Leakage(Ntechnologies+1,k) = sum(leakList);
    % Add new leaks to the list of generated leaks
    leakList = [leakList; newLeaks(k).Q];
    %% Tech modules: loop through each LDAR program
    for i = 1: Ntechnologies
        Leakage(i,k) = sum(leakStruct(i).Q);
        %Define a function handle for each technology's detection function
        fh = str2func([legendString{i} 'Detect']);
        [leaksFound(i).Q,leakStruct(i),repairCost(i,k),financeStruct{i}] = fh(techStruct{i},financeStruct{i},processStruct{i},leakStruct(i),leaksFound(i).Q,repairdata,timeStruct,k,atmosphere);
        % Remove leaks through the null process
        [nullRepairedIndeces, remainingIndeces] = NullLeakFixer(length(leakStruct(i).Q),GasField.NullRepairRate,deltaT);
        nullRepaired(i).Q = [nullRepaired(i).Q; leakStruct(i).Q(nullRepairedIndeces)];
        leakStruct(i) = structfun(@(M)M(remainingIndeces),leakStruct(i),'Uniform',0);
        repairCost(i,k) = repairCost(i,k) + sum(datasample(repairdata,length(nullRepairedIndeces)));
        % Now add the new leaks
        fields = fieldnames(leakStruct(i));
        for j = 1:numel(fields)
            leakStruct(i).(fields{j}) = [leakStruct(i).(fields{j}); newLeaks(k).(fields{j})];
        end
    end
end

%% Save results
clear('i','k','j','atmosphere','endTime','deltaT','fh','fields','nullRepairedIndeces','remainingIndeces','repairdata','time','timeSteps')
[Capital, Maintenance, Finding] = FinanceSummary(financeStruct);

if ~exist('Results','dir')
mkdir([pwd '/Results'])
end
save([pwd '/Results/Realization' num2str(length(dir([pwd '/Results']))-1) '.mat'])
