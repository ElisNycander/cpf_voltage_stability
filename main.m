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
CPFOptions.windGenerators = [1]; % generators where to put wind farms 
CPFOptions.windBusType = 'pv'; % type of wind production: pq or pv

CPFOptions.powerFactor = 0.8;
CPFOptions.pWind = 200:100:200; % wind capacity

%% runcpfs options - contingencies and target case
CPFOptions.caseFile = 'case4gs';

CPFOptions.tripAllLines = 0; % trip all lines; overrides tripLines
CPFOptions.tripLines = [1 2]; % list of lines to trip
CPFOptions.tripGenerators = []; % list of generators to trip

CPFOptions.loadIncreaseBuses = [3]; % if empty the load is increased at all buses
CPFOptions.productionIncreaseGenerators = []; % if empty all load is compensated at slack bus

%% plotCPF OPTIONS
CPFOptions.Plot = struct();

%% stability margin plots
CPFOptions.Plot.stabilityMargin = 1; % plot stability margin curves
CPFOptions.Plot.separateContingencies = 0; % separate Plot for each contingency

%% plots for invidiual cases
CPFOptions.Plot.individualCases = 1; % plot voltages and power curves for invidual scenarios

CPFOptions.Plot.pvCurve = 1; % Plot pvCurve curves
CPFOptions.Plot.pqCurve = 1; % Plot pqCurve curves
CPFOptions.Plot.ppCurve = 1; % Plot ppCurve curves

CPFOptions.Plot.voltageBuses = [1 2 3 4]; % buses to Plot voltage for
CPFOptions.Plot.powerBuses = [1 2 3 4]; % buses to Plot Q for
CPFOptions.Plot.powerInjectedOrGenerated = 'gen'; % gen: Plot Q generation, inj: Plot Q injection (includes load)
CPFOptions.Plot.subfigures = 1; % all plots in same figure

CPFOptions.Plot.colormap = 'lines';




%% options for cpf
mpopt = mpoption();
mpopt.verbose = 0;
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

