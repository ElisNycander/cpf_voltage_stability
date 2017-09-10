

function cpfScenarios = runCPF(CPFOptions,mpopt)

define_constants;

% base case
%mpcb = CPFOptions.caseFunction();
mpcb = CPFOptions.mpc;

% add WF to gen matrix
if strcmp(CPFOptions.windScheme,'buses')
   % add wind generator for each bus
   grow = size(mpcb.gen,1);
   gcol = size(mpcb.gen,2);
   nwf = length(CPFOptions.windBuses);
   mpcb.gen = [mpcb.gen ; zeros(nwf,gcol)];
   mpcb.gen(grow+1:grow+nwf,GEN_BUS) = CPFOptions.windBuses; % fill in generator bus
   
   % also add rows to gencost matrix
   if isfield(mpcb,'gencost')
       mpcb.gencost = [mpcb.gencost ; zeros(nwf,size(mpcb.gencost,2))];
   end
       
   for i=1:nwf
       % set WF gen voltage to conform with other generators at the same bus
       gi = find(mpcb.gen(1:grow,GEN_BUS)==CPFOptions.windBuses(i));
       
       mpcb.gen(grow+i,MBASE) = mpcb.baseMVA;
       if ~isempty(gi)
           gi = gi(1);
           mpcb.gen(grow+i,[VG]) = mpcb.gen(gi,[VG]); % WF gen voltage = 1
       else
           mpcb.gen(grow+i,VG) = 1;
       end
   end
   
   % for WF pv type, set all WF buses to pv and activate WF generators
   if strcmp(CPFOptions.windBusType,'pv')
       mpcb.bus(CPFOptions.windBuses,BUS_TYPE) = 2;
       mpcb.gen(grow+1:grow+nwf,GEN_STATUS) = 1;
   end
   
   if CPFOptions.removeOtherGeneration
       for i=1:length(CPFOptions.windBuses)
            % deactive old generators at WF buses
            gentmp = mpcb.gen(1:grow,:);
            gentmp( gentmp(:,GEN_BUS) == CPFOptions.windBuses(i) , GEN_STATUS) = 0;
            mpcb.gen(1:grow,:) = gentmp;
       end  
       if strcmp(CPFOptions.windBusType,'pq')
            % buses that were PV and had generators deactivated should now
            % be PQ
            mpcb.bus( CPFOptions.windBuses,BUS_TYPE) = 1;
       end
   end
end

% save base case for later
mpcbSaved = mpcb;

cpfScenarios = cell(CPFOptions.nContingencies,CPFOptions.nWindPoints); % for storing results of cpfs

for k=1:CPFOptions.nWindPoints
    
    mpcb = mpcbSaved;
    
    %% add wind power [P Q] to gen matrix
    % If wind farms are PQ buses [PG QMAX] is only used as storage for
    % negative load. If wind farms are PV buses [PG QMAX] are used for
    % corresponding PV bus.
    
    if strcmp(CPFOptions.windScheme,'buses')
        
        if length(CPFOptions.windBusShare) ~= length(nwf)
            windWeights = 1/length(nwf) * ones(1,nwf); % proportional share at all buses
        else
            windWeights = CPFOptions.windBusShare / sum(CPFOptions.windBusShare);
        end
        
        for iii=1:nwf
            
            qfactor = sqrt(1-CPFOptions.powerFactor^2)/CPFOptions.powerFactor;
            if strcmp(CPFOptions.powerAngle,'lead')
                qfactor = -qfactor;
            end
            mpcb.gen(iii+grow,[PG QMAX QMIN]) = [CPFOptions.pWind(k) ... % set [PG QMAX QMIN]
                CPFOptions.pWind(k)*qfactor ...
                -CPFOptions.pWind(k)*qfactor ] * windWeights(iii);
        end
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
        for i=1:nwf
            % find bus of generator
            bi = CPFOptions.windBuses(i);
            gi = grow + i;
            % same constant power factor load in base and target case
            mpcb.bus(bi,[PD QD]) = mpcb.bus(bi,[PD QD]) - mpcb.gen(gi,[PG QMAX]);
            mpct.bus(bi,[PD QD]) = mpct.bus(bi,[PD QD]) - mpcb.gen(gi,[PG QMAX]);
            
        end
    end
    
    
    
    for i=1:CPFOptions.nContingencies
        
        cpfScenarios{i,k} = struct();
        mpc1 = mpcb;
        mpc2 = mpct;
        
        %% apply contingency
        if i == 1 % no fault
            cpfScenarios{i,k}.str = sprintf(['No fault\n%i MW wind (' CPFOptions.windBusType ')'],CPFOptions.pWind(k));
        elseif i <= CPFOptions.nLines+1 % trip line
            idx = CPFOptions.tripLines(i-1);
            mpc1.branch(idx,:) = [];
            mpc2.branch(idx,:) = [];
            cpfScenarios{i,k}.str = sprintf(['Trip line from bus %i to bus %i\n%i MW wind (' CPFOptions.windBusType ')'], ...
                [mpcb.branch(idx,[1 2]) CPFOptions.pWind(k)]);
            
        elseif i <= CPFOptions.nLines+1 + CPFOptions.nGenerators  % trip generator
            idx = CPFOptions.tripGenerators(i-CPFOptions.nLines-1);
            mpc1.gen(idx,GEN_STATUS) = 0;
            mpc2.gen(idx,GEN_STATUS) = 0;
            cpfScenarios{i}.str = sprintf(['Trip generator at bus %i\n%i MW wind (' CPFOptions.windBusType ')'], ...
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
            
            
            if CPFOptions.mergeGeneration
                mpcm1 = mergeGenerators(mpc1);
                mpcm2 = mergeGenerators(mpc2);
                [mpc3,success] = runcpf(mpcm1,mpcm2,mpopt);
            else
                [mpc3,success] = runcpf(mpc1,mpc2,mpopt);
            end
            iter = iter + 1;
        end
        %mpc3.cpf.loadIncreaseIterations = iter; % needed to know original base load
        
        if strcmp(CPFOptions.windBusType,'pq')
            % remove negative load generation
            for ii=1:nwf
                % find bus of generator
                bi = CPFOptions.windBuses(i);
                gi = grow + i;
                
                mpc1.bus(bi,[PD QD]) = mpc1.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QMAX]);
                mpc2.bus(bi,[PD QD]) = mpc2.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QMAX]);
                mpc3.bus(bi,[PD QD]) = mpc3.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QMAX]);
            end
        end
        
        cpfScenarios{i,k}.mpcb = mpc1;
        cpfScenarios{i,k}.mpct = mpc2;
        cpfScenarios{i,k}.mpcc = mpc3;
    end
end