

function cpfScenarios = runCPF(CPFOptions,mpopt)

define_constants;

% base case
mpcb = CPFOptions.caseFunction();

if strcmp(CPFOptions.windBusType,'pq') % change wind farm buses to PQ
    for i=1:length(CPFOptions.windGenerators)
        gi = CPFOptions.windGenerators(i);
        bi = mpcb.gen(gi, GEN_BUS);
        
        mpcb.gen(gi, GEN_STATUS) = 0; % deactivate generator
        mpcb.bus(bi, BUS_TYPE) = 1; % change to PQ node
    end
end

% save base case for later
mpcbSaved = mpcb;

cpfScenarios = cell(CPFOptions.nContingencies,CPFOptions.nWindPoints); % for storing results of cpfs

for k=1:CPFOptions.nWindPoints
    
    mpcb = mpcbSaved;
    
    %% add wind power generation
    % If wind farms are PQ buses [PG QMAX] is only used as storage for
    % negative load. If wind farms are PV buses [PG QMAX] are used for
    % corresponding PV bus.
    
    for iii=1:length(CPFOptions.windGenerators)
        qfactor = sqrt(1-CPFOptions.powerFactor^2)/CPFOptions.powerFactor;
        if strcmp(CPFOptions.powerAngle,'lead')
            qfactor = -qfactor;
        end
        mpcb.gen(iii,[PG QMAX QMIN]) = [CPFOptions.pWind(k) ... % change [PG QMAX]
                                        CPFOptions.pWind(k)*qfactor ...
                                       -CPFOptions.pWind(k)*qfactor ];
    end
    
    
    %% form target load case
    mpct = mpcb;
    
    % proportional load increase at specified buses
    activeLoadIncreaseVector = mpcb.bus(CPFOptions.loadIncreaseBuses,PD);
    reactiveLoadIncreaseVector = mpcb.bus(CPFOptions.loadIncreaseBuses,QD);
    activeLoadIncrease = sum(activeLoadIncreaseVector);
    
    if CPFOptions.loadIterations && CPFOptions.fromZeroLoad % begin first CPF from 0 load
        mpcb.bus(CPFOptions.loadIncreaseBuses,[PD QD]) =  ...
            mpcb.bus(CPFOptions.loadIncreaseBuses,[PD QD]) - ...
            [activeLoadIncreaseVector reactiveLoadIncreaseVector];
    else % begin CPF from base case load
        mpct.bus(CPFOptions.loadIncreaseBuses,[PD QD]) =  ...
            mpct.bus(CPFOptions.loadIncreaseBuses,[PD QD]) + ...
            [activeLoadIncreaseVector reactiveLoadIncreaseVector];
    end
    
    % QUESTION: CAN QLIM CHANGE GRADUALLY IN CPF?? - NO
    % increase generation
    if ~isempty(CPFOptions.productionIncreaseGenerators)
        % add up generation of specified generators
        gi = CPFOptions.productionIncreaseGenerators;
        totalGeneration = sum(mpcb.gen(gi,PG));
        scaleFactor = activeLoadIncrease / totalGeneration;
        % Note: For PV nodes, QMAX is constant in CPF (value from base case scenario?)
        % so only the active power set point is changed.
        mpct.gen(gi,[PG QMAX]) = mpct.gen(gi,[PG QMAX]) * scaleFactor;
    end
    
    if strcmp(CPFOptions.windBusType,'pq') % add negative load generation
        for i=1:length(CPFOptions.windGenerators)
            % find bus of generator
            gi = CPFOptions.windGenerators(i);
            bi = mpcb.gen(1,gi);
            
            mpcb.bus(bi,[PD QD]) = mpcb.bus(bi,[PD QD]) - mpcb.gen(gi,[PG QMAX]);
            mpct.bus(bi,[PD QD]) = mpct.bus(bi,[PD QD]) - mpct.gen(gi,[PG QMAX]);
            
        end
    end
    
    
    for i=1:CPFOptions.nContingencies
        
        cpfScenarios{i,k} = struct();
        mpc1 = mpcb;
        mpc2 = mpct;
        
        %% apply contingency
        if i == 1 % no fault
            cpfScenarios{i,k}.str = sprintf(['No fault\n%i MW wind'],CPFOptions.pWind(k));
        elseif i <= CPFOptions.nLines+1 % trip line
            idx = CPFOptions.tripLines(i-1);
            mpc1.branch(idx,:) = [];
            mpc2.branch(idx,:) = [];
            cpfScenarios{i,k}.str = sprintf(['Trip line from bus %i to bus %i\n%i MW wind'], ...
                [mpcb.branch(idx,[1 2]) CPFOptions.pWind(k)]);
            
        elseif i <= CPFOptions.nLines+1 + CPFOptions.nGenerators  % trip generator
            idx = CPFOptions.tripGenerators(i-CPFOptions.nLines-1);
            mpc1.gen(idx,GEN_STATUS) = 0;
            mpc2.gen(idx,GEN_STATUS) = 0;
            cpfScenarios{i}.str = sprintf('Trip generator at bus %i\n%i MW wind', ...
                [mpcb.gen(idx,1) CPFOptions.pWind(k)]);
        end
        
        iter = 0;
        success = 0;
        baseLoad = mpc1.bus(CPFOptions.loadIncreaseBuses,[PD QD]);
        targetLoad = mpc2.bus(CPFOptions.loadIncreaseBuses,[PD QD]);
        while ~success && iter < CPFOptions.maxLoadIterations 
            % increase base case load
            mpc1.bus(CPFOptions.loadIncreaseBuses,[PD QD]) = baseLoad + ...
                [activeLoadIncreaseVector reactiveLoadIncreaseVector] * iter;
            
            mpc2.bus(CPFOptions.loadIncreaseBuses,[PD QD]) = targetLoad + ...
                [activeLoadIncreaseVector reactiveLoadIncreaseVector] * iter;
            
            [mpc3,success] = runcpf(mpc1,mpc2,mpopt);
            iter = iter + 1;
        end
        %mpc3.cpf.loadIncreaseIterations = iter; % needed to know original base load
        
        if strcmp(CPFOptions.windBusType,'pq')
            % remove negative load generation
            for ii=1:length(CPFOptions.windGenerators)
                % find bus of generator
                gi = CPFOptions.windGenerators(ii);
                bi = mpcb.gen(1,gi);
                
                mpc1.bus(bi,[PD QD]) = mpc1.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QG]);
                mpc2.bus(bi,[PD QD]) = mpc2.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QG]);
                mpc3.bus(bi,[PD QD]) = mpc3.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QG]);
            end
        end
        
        cpfScenarios{i,k}.mpcb = mpc1;
        cpfScenarios{i,k}.mpct = mpc2;
        cpfScenarios{i,k}.mpcc = mpc3;
    end
end