%% find voltage stability limits
%
%
%
%
%
%%

clear;
close all;
define_constants;
CPFOptions = struct();

%% general options
CPFOptions.filename = 'pv_run';

%% wind power options

% windScheme
% buses - Specify buses with WF and total amount and proportional allocation.
% WFs may be PQ (negative load) or PV. 
% generators - Specify generators where to increase WF. WFs may be PQ or PV. 
CPFOptions.windScheme = 'buses';
CPFOptions.removeOtherGeneration = 0; % 1 - remove other generation at WF buses
CPFOptions.windGenerators = []; % generators where to put wind farms 
CPFOptions.windBuses = [4]; % buses where to put WFs
CPFOptions.windBusShare = []; % ratio of wind farm share
CPFOptions.windBusType = 'pq'; % type of wind production: pq or pv

powerAngle = 35.53; % 30 35.53 % 40.26 45
CPFOptions.powerFactor = cos(powerAngle/180*pi);
%CPFOptions.powerFactor = 0.9;
CPFOptions.powerAngle = 'lag'; % lead/lag 
CPFOptions.pWind = 1900; % wind capacity

%% runcpfs options - contingencies and target case
CPFOptions.caseFile = 'case4gs';

% The maximum loadability will be calculated as the fist point where 
% voltages start increasing, as the CPF solver may take a few false steps at 
% the end with icreasing voltages, presumably due to numerical issues
CPFOptions.checkCpfTermination = 1; 
CPFOptions.voltageTolerance = 1e-3;

CPFOptions.loadIterations = 1; % option to run CPF with different base load if it does not converge
CPFOptions.fromZeroLoad = 0; % first CPF from zero load, otherwise begin at load given in case  
CPFOptions.maxLoadIterations = 5; % maximum number of attempts to increase load

CPFOptions.tripAllLines = 0; % trip all lines; overrides tripLines
CPFOptions.tripLines = []; % list of lines to trip
CPFOptions.tripGenerators = []; % list of generators to trip

CPFOptions.loadIncreaseBuses = [3]; % if empty the load is increased at all buses
CPFOptions.productionIncreaseGenerators = []; % if empty all load is compensated at slack bus

%% plotCPF OPTIONS
CPFOptions.Plot = struct();

%% stability margin plots
CPFOptions.Plot.stabilityMargin = 1; % plot stability margin curves
CPFOptions.Plot.separateContingencies = 1; % separate Plot for each contingency

%% plots for invidiual cases
CPFOptions.Plot.individualCases = 1; % plot voltages and power curves for invidual scenarios

CPFOptions.Plot.pvCurve = 1; % Plot pvCurve curves
CPFOptions.Plot.pqCurve = 1; % Plot pqCurve curves
CPFOptions.Plot.ppCurve = 1; % Plot ppCurve curves

CPFOptions.Plot.voltageBuses = [1 2 3 4]; % buses to Plot voltage for
CPFOptions.Plot.powerBuses = [1 2 3 4]; % buses to Plot Q for
CPFOptions.Plot.powerInjectedOrGenerated = 'gen'; % gen: Plot Q generation, inj: Plot Q injection (includes load)
CPFOptions.Plot.subfigures = 0; % all plots in same figure

CPFOptions.Plot.colormap = 'lines';




%% options for cpf
mpopt = mpoption();
mpopt.verbose = 1;
mpopt.cpf.enforce_q_lims = 1;
mpopt.out.all = 0;
mpopt.cpf.user_callback = 'cpf_modified_callback';
mpopt.cpf.parameterization = 3;
mpopt.adapt_step = 0;
mpopt.step = 0.1;


%% do some calculations
CPFOptions = processCPFOptions(CPFOptions);

%% run cpfs

% run cpfs
tic;      
CPFResults = runCPF(CPFOptions,mpopt);
display(fprintf('Ran %i CPFs in %0.2f seconds',[size(CPFResults,1)*size(CPFResults,2), toc]));

% collect results into matrices
[processedCPFResults,CPFResults] = processCPF(CPFResults,CPFOptions);
    
% plot results
plotCPF(processedCPFResults,CPFResults,CPFOptions);

% save results
results = struct('processedCPFResults',processedCPFResults, ...
                 'CPFResults',CPFResults, ...
                 'CPFOptions',CPFOptions );
save(CPFOptions.filename,'results');

