

function contingencyCases = runCPF(mpcb,CPFOptions,mpopt)

define_constants;

%% form target load case
mpct = mpcb;

% proportional load increase at specified buses
mpct.bus(CPFOptions.loadIncreaseBuses,[PD QD]) =  ...
    mpct.bus(CPFOptions.loadIncreaseBuses,[PD QD]) * 2;
loadIncrease = sum(mpcb.bus(CPFOptions.loadIncreaseBuses,PD));



% QUESTION: CAN QLIM CHANGE GRADUALLY IN CPF?? - NO
% increase generation
if ~isempty(CPFOptions.productionIncreaseGenerators)
    % add up generation of specified generators
    gi = CPFOptions.productionIncreaseGenerators;
    totalGeneration = sum(mpcb.gen(gi,PG));
    scaleFactor = loadIncrease / totalGeneration;
    % Note: For PV nodes, QMAX is constant in CPF (value from base case scenario?)
    % so only the active power set point is changed. 
    mpct.gen(gi,[PG QMAX]) = mpct.gen(gi,[PG QMAX]) * scaleFactor;
end

if strcmp(CPFOptions.windBusType,'pq')
    % add negative load generation
    for i=1:length(CPFOptions.windGenerators)
        % find bus of generator
        gi = CPFOptions.windGenerators(i);
        bi = mpcb.gen(1,gi);
        
        mpcb.bus(bi,[PD QD]) = mpcb.bus(bi,[PD QD]) - mpcb.gen(gi,[PG QMAX]);
        mpct.bus(bi,[PD QD]) = mpct.bus(bi,[PD QD]) - mpct.gen(gi,[PG QMAX]);
        
    end
end

contingencyCases = cell(1,CPFOptions.nContingencies); % for storing results of cpfs

for i=1:CPFOptions.nContingencies
    
    contingencyCases{i} = struct();
    mpc1 = mpcb;
    mpc2 = mpct;
    
    %% apply contingency
    if i == 1 % no fault
        contingencyCases{i}.str = 'No fault';
    elseif i <= CPFOptions.nLines+1 % trip line
        idx = CPFOptions.tripLines(i-1);
        mpc1.branch(idx,:) = [];
        mpc2.branch(idx,:) = [];
        contingencyCases{i}.str = sprintf('Trip line from bus %i to bus %i',mpcb.branch(idx,[1 2]));
    elseif i <= CPFOptions.nLines+1 + CPFOptions.nGenerators  % trip generator
        idx = CPFOptions.tripGenerators(i-CPFOptions.nLines-1);
        mpc1.gen(idx,GEN_STATUS) = 0;
        mpc2.gen(idx,GEN_STATUS) = 0; 
        contingencyCases{i}.str = sprintf('Trip generator at bus %i',mpcb.gen(idx,1));
    end
    
    
    [mpc3,~] = runcpf(mpc1,mpc2,mpopt);
    
    if strcmp(CPFOptions.windBusType,'pq')
        % remove negative load generation
        for ii=1:length(CPFOptions.gen_negative_load)
            % find bus of generator
            gi = CPFOptions.gen_negative_load(ii);
            bi = mpcb.gen(1,gi);
            
            mpc1.bus(bi,[PD QD]) = mpc1.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QG]);
            mpc2.bus(bi,[PD QD]) = mpc2.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QG]);
            mpc3.bus(bi,[PD QD]) = mpc3.bus(bi,[PD QD]) + mpcb.gen(gi,[PG QG]);
        end
    end

    contingencyCases{i}.mpcb = mpc1;
    contingencyCases{i}.mpct = mpc2;
    contingencyCases{i}.mpcc = mpc3;
end