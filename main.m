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

%% wind power options
CPFOptions.windGenerators = [1]; % generators where to put wind farms 
CPFOptions.windBusType = 'pv'; % type of wind production: pq or pv

CPFOptions.powerFactor = 0.9;
CPFOptions.pWind = 0:100:2000; % wind capacity

%% runcpfs options - contingencies and target case
CPFOptions.caseFile = 'case4gs';

CPFOptions.tripAllLines = 1; % trip all lines; overrides tripLines
CPFOptions.tripLines = []; % list of lines to trip
CPFOptions.tripGenerators = []; % list of generators to trip

CPFOptions.loadIncreaseBuses = [3]; % if empty the load is increased at all buses
CPFOptions.productionIncreaseGenerators = []; % if empty all load is compensated at slack bus


%% plotcpfs options
CPFOptions.Plot = struct();

CPFOptions.Plot.disable = 1;

CPFOptions.Plot.pvCurve = 1; % Plot pvCurve curves
CPFOptions.Plot.pqCurve = 1; % Plot pqCurve curves
CPFOptions.Plot.ppCurve = 1; % Plot ppCurve curves
CPFOptions.Plot.voltageStabilityBar = 1; % Plot voltage stability bar Plot

CPFOptions.Plot.voltageBuses = [1 2 3 4]; % buses to Plot voltage for
CPFOptions.Plot.powerBuses = [1 4]; % buses to Plot Q for
CPFOptions.Plot.powerInjectedOrGenerated = 'gen'; % gen: Plot Q generation, inj: Plot Q injection (includes load)
CPFOptions.Plot.subfigures = 1; % all plots in same figure

CPFOptions.Plot.colormap = 'lines';

CPFOptions.Plot.separateContingencies = 0; % separate Plot for each contingency (stability margin voltageStabilityBar CPFOptions.pWind)


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

% base case
mpcb = CPFOptions.caseFunction();

if strcmp(CPFOptions.windBusType,'pqCurve')
    % change wind farm buses to PQ
    for i=1:length(CPFOptions.windGenerators)
        gi = CPFOptions.windGenerators(i);
        bi = mpcb.gen(gi, GEN_BUS);
        
        mpcb.gen(gi, GEN_STATUS) = 0; % deactivate generator
        mpcb.bus(bi, BUS_TYPE) = 1; % change to PQ node
    end
end

N = CPFOptions.nWindPoints;
powerFactor = CPFOptions.powerFactor;

Psecure = zeros(CPFOptions.nContingencies,N);
Pmax = zeros(CPFOptions.nContingencies,N);
secure = zeros(CPFOptions.nContingencies,N);
pMaxNminus1 = zeros(1,N);
pSecureNminus1 = zeros(1,N);
Nminus1 = zeros(1,N);


tic;
for ii=1:N % increasinge wind capacity
    
    
    % If wind farms are PQ buses [PG QMAX] is only used as storage for
    % negative load. If wind farms are PV buses [PG QMAX] are used for
    % corresponding PV bus. 
    for iii=1:length(CPFOptions.windGenerators)
        mpcb.gen(iii,[PG QMAX QMIN]) = [CPFOptions.pWind(ii) ... % change [PG QMAX]
            CPFOptions.pWind(ii)*sqrt(1-powerFactor^2)/powerFactor ...
            CPFOptions.pWind(ii)*sqrt(1-powerFactor^2)/powerFactor ];
    end
    
    c = runCPF(mpcb,CPFOptions,mpopt);

    [cpf_res,c] = processCPF(c,CPFOptions);
    
    Psecure(:,ii) = cpf_res.Psecure';
    Pmax(:,ii) = cpf_res.Pmax';
    secure(:,ii) = cpf_res.secure';
    pMaxNminus1(ii) = cpf_res.pMaxNminus1;
    pSecureNminus1(ii) = cpf_res.pSecureNminus1;
    Nminus1(ii) = cpf_res.Nminus1;
end
display(sprintf('Ran %i CPFs in %0.2f seconds',[size(Pmax,1)*size(Pmax,2), toc]));
contingencies = cpf_res.contingencies;


CPFResults = struct( ...
                    'pMax',Pmax, ...
                    'pSecure',Psecure, ...
                    'pMaxNminus1',pMaxNminus1, ...
                    'pSecureNminus1',pSecureNminus1, ...
                    'Nminus1',Nminus1, ...
                    'pWind',CPFOptions.pWind ...
);
CPFResults.contingencies = contingencies;

plotCPF(CPFResults,CPFOptions);

