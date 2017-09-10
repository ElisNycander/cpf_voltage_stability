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
CPFOptions.removeOtherGeneration = 1; % 1 - remove other generation at WF buses

CPFOptions.windBuses = [5]; % buses where to put WFs
CPFOptions.windBusShare = []; % ratio of wind farm share
CPFOptions.windGenerators = []; % generators where to put wind farms 
CPFOptions.windBusType = 'pv'; % type of wind production: pq or pv

qpRatio = 0.2;
powerAngle=atan(qpRatio);
CPFOptions.powerFactor = cos(powerAngle/180*pi);
CPFOptions.powerAngle = 'lag'; % lead/lag 
CPFOptions.pWind = 450; % wind capacity

%% runcpfs options - contingencies and target case
CPFOptions.caseFile = 'case24_ieee_rts';

CPFOptions.mergeGeneration = 1;

% The maximum loadability will be calculated as the fist point where 
% voltages start increasing, as the CPF solver may take a few false steps at 
% the end with icreasing voltages, presumably due to numerical issues
CPFOptions.checkCpfTermination = 0; 
CPFOptions.voltageTolerance = 1e-3;

CPFOptions.loadIterations = 0; % option to run CPF with different base load if it does not converge
CPFOptions.fromZeroLoad = 0; % first CPF from zero load, otherwise begin at load given in case  
CPFOptions.maxLoadIterations = 3; % maximum number of attempts to increase load

CPFOptions.tripAllLines = 0; % trip all lines; overrides tripLines
CPFOptions.tripLines = []; % list of lines to trip
CPFOptions.tripGenerators = []; % list of generators to trip

CPFOptions.loadIncreaseBuses = [11:24]; % if empty the load is increased at all buses
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

CPFOptions.Plot.voltageBuses = []; % buses to Plot voltage for
CPFOptions.Plot.powerBuses = [1 2 7 13 15 16 18 21 22 23]; % buses to Plot Q for
%CPFOptions.Plot.powerBuses = [13];
CPFOptions.Plot.powerInjectedOrGenerated = 'gen'; % gen: Plot Q generation, inj: Plot Q injection (includes load)
CPFOptions.Plot.subfigures = 0; % all plots in same figure

CPFOptions.Plot.colormap = 'lines';




%% options for cpf
mpopt = mpoption();
mpopt.verbose = 2;
mpopt.cpf.enforce_q_lims = 1;
mpopt.out.all = 0;
mpopt.cpf.user_callback = 'cpf_modified_callback';
mpopt.cpf.parameterization = 3;
mpopt.adapt_step = 0;
mpopt.step = 0.1;


%% do some calculations
CPFOptions = processCPFOptions(CPFOptions);

%% set Q limits

CPFOptions.mpc.gen(:,[QMIN]) = [-1e9];
CPFOptions.mpc.gen(CPFOptions.mpc.gen(:,GEN_BUS) == 13,QMAX) = 1e9;
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

