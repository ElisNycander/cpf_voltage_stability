clear;
close all;
define_constants;
CPFOptions = struct();

%% script options
% options particular to this script


%% wind power options

% windScheme
% buses - Specify buses with WF and total amount and proportional allocation.
% WFs may be PQ (negative load) or PV. 
% generators - Specify generators where to increase WF. WFs may be PQ or PV. 
CPFOptions.windScheme = 'buses';
CPFOptions.removeOtherGeneration = 1; % 1 - remove other generation at WF buses
CPFOptions.windGenerators = []; % generators where to put wind farms 
CPFOptions.windBuses = [4]; % buses where to put WFs
CPFOptions.windBusShare = []; % ratio of wind farm share
CPFOptions.windBusType = 'pv'; % type of wind production: pq or pv

CPFOptions.powerFactor = 0.8;
CPFOptions.powerAngle = 'lag'; % lead/lag 
CPFOptions.pWind = 0:100:2000; % wind capacity

%powerAngle = linspace(0,pi/4,20);
pqRatio = linspace(0,1,21);
powerAngle = atan(pqRatio);
%powerAngle = linspace(-pi/8,pi/8,2);
%powerAngle = pi/180*[30 35];
powerFactor = cos(powerAngle);

%% runcpfs options - contingencies and target case
CPFOptions.caseFile = 'case4gs';

CPFOptions.mergeGeneration = 1; % use only one generator per bus in cpf

% The maximum loadability will be calculated as the fist point where 
% voltages start increasing, as the CPF solver may take a few false steps at 
% the end with icreasing voltages, presumably due to numerical issues
CPFOptions.checkCpfTermination = 0; 
CPFOptions.voltageTolerance = 1e-3;

CPFOptions.loadIterations = 1; % option to run CPF with different base load if it does not converge
CPFOptions.fromZeroLoad = 0; % first CPF from zero load, otherwise begin at load given in case  
CPFOptions.maxLoadIterations = 5; % maximum number of attempts to increase load

CPFOptions.tripAllLines = 0; % trip all lines; overrides tripLines
CPFOptions.tripLines = []; % list of lines to trip
CPFOptions.tripGenerators = []; % list of generators to trip

CPFOptions.loadIncreaseBuses = []; % if empty the load is increased at all buses
CPFOptions.productionIncreaseGenerators = []; % if empty all load is compensated at slack bus

%% plotCPF OPTIONS
CPFOptions.Plot = struct();

%% stability margin plots
CPFOptions.Plot.stabilityMargin = 0; % plot stability margin curves
CPFOptions.Plot.separateContingencies = 1; % separate Plot for each contingency

%% plots for invidiual cases
CPFOptions.Plot.individualCases = 0; % plot voltages and power curves for invidual scenarios

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

    
    
    
    % do some calculations
    CPFOptions = processCPFOptions(CPFOptions);
    
    
    CPFOptions.mpc.gen(2,[QMAX QMIN]) = [318 -1e9];
    CPFOptions.mpc.gen(1,[QMIN]) = [-1e9];
    
    nPowerFactors = length(powerFactor);
    
    loadabilityMargin = zeros(nPowerFactors,CPFOptions.nWindPoints);
    securityMargin = zeros(nPowerFactors,CPFOptions.nWindPoints);
    securityLimitType = zeros(nPowerFactors,CPFOptions.nWindPoints);
    basePLossFraction = zeros(size(securityMargin));
    nosePLossFraction = zeros(size(securityMargin));
     
for ii=1:nPowerFactors
    
    CPFOptions.powerFactor = powerFactor(ii);
    if powerAngle(ii) > 0
        CPFOptions.powerAngle = 'lag';
    else
        CPFOptions.powerAngle = 'lead';
    end
    % run cpfs
    tic;
    CPFResults = runCPF(CPFOptions,mpopt);
    display(sprintf('Ran %i CPFs in %0.2f seconds',[size(CPFResults,1)*size(CPFResults,2), toc]));
    
    % collect results into matrices
    [results,] = processCPF(CPFResults,CPFOptions);
    
    % plot stability margins
    plotCPF(results,struct(),CPFOptions);
    
    loadabilityMargin(ii,:) = results.pMax;
    securityMargin(ii,:) = results.pSecure;
    securityLimitType(ii,:) = results.securityLimitType;
    
    basePLossFraction(ii,:) = results.basePLossFraction;
    nosePLossFraction(ii,:) = results.nosePLossFraction;
end

results  = struct('powerFactor',powerFactor, ...
                    'pWind',CPFOptions.pWind, ...
                    'loadabilityMargin',loadabilityMargin, ...
                    'securityMargin',securityMargin, ...
                    'securityLimitType',securityLimitType, ...
                    'powerAngle',powerAngle, ...
                    'basePLossFraction',basePLossFraction, ...
                    'nosePLossFraction',nosePLossFraction);

save('pf_results.mat','results');

plotPFSurface('pf_results.mat');

