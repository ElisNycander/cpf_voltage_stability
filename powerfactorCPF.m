clear;
close all;
define_constants;
CPFOptions = struct();

%% wind power options
CPFOptions.windGenerators = [1]; % generators where to put wind farms 
CPFOptions.windBusType = 'pq'; % type of wind production: pq or pv

CPFOptions.powerFactor = 0.8;
CPFOptions.powerAngle = 'lag'; % lead/lag 
CPFOptions.pWind = 0:100:2000; % wind capacity

powerAngle = linspace(-pi/4,pi/4,20);
powerFactor = cos(powerAngle);

%% runcpfs options - contingencies and target case
CPFOptions.caseFile = 'case4gs';

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


%% do some calculations
CPFOptions = processCPFOptions(CPFOptions);

nPowerFactors = length(powerFactor);

loadabilityMargin = zeros(nPowerFactors,CPFOptions.nWindPoints);
securityMargin = zeros(nPowerFactors,CPFOptions.nWindPoints);

for i=1:nPowerFactors
    
    CPFOptions.powerFactor = powerFactor(i);
    if powerAngle(i) > 0
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
    
    loadabilityMargin(i,:) = results.pMax;
    securityMargin(i,:) = results.pSecure;
end

results  = struct('powerFactor',powerFactor, ...
                    'pWind',CPFOptions.pWind, ...
                    'loadabilityMargin',loadabilityMargin, ...
                    'securityMargin',securityMargin );

figure;
surf(results.pWind,180/pi*powerAngle,loadabilityMargin);
title('Loadability margin')
xlabel('Wind power (MW)');
ylabel('Power angle (degrees)');
zlabel('Loadability margin (MW)');
figure;
surf(results.pWind,180/pi*powerAngle,securityMargin);
title('Security margin');
xlabel('Wind power (MW)');
ylabel('Power angle (degrees)');
zlabel('Security margin (MW)');


