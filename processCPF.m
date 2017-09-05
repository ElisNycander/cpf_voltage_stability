% [res,contingencyCases] = processCPF(contingencyCases,CPFOptions)
% Collect and Plot results from cpfs contained in cell array contingencyCases
%
% The following extra results are stored in each mpcc.cpf:
% Sinj(bus x step) - net injected complex power at each bus
% Pgen(bus x step) - generated active power (injected power + load)
% Qgen(bus x step) - generated reactive power (injected power + load)
% pMax             - system total load at maximum loadability point (0 if
%                    base case was not solvable)
% pSecure          - system total load at last point satisfying operating limits
%

function [res,contingencyCases] = processCPF(contingencyCases,CPFOptions)

define_constants;

nContingencies = size(contingencyCases,1);
nWindPoints = size(contingencyCases,2);

pMax = zeros(nContingencies,nWindPoints);
pSecure = zeros(nContingencies,nWindPoints);
max_lam = zeros(nContingencies,nWindPoints); % maximum lamda when stability limit is reached (system unsolvable)
secure_lam = zeros(nContingencies,nWindPoints); % maximum lambda for which system is secure
contingencies = cell(nContingencies,1);
secure = zeros(nContingencies,nWindPoints); % binary = 1 if base case power flow secure
securityLimitType = zeros(nContingencies,nWindPoints); % type of security violation, 1 min, 2 max, 0 no violation, -1 no convergence

for k=1:nWindPoints % loop over wind power
    for i=1:nContingencies % loop contingencies
        
        mpcb = contingencyCases{i,k}.mpcb;
        mpct = contingencyCases{i,k}.mpct;
        mpcc = contingencyCases{i,k}.mpcc;
        
        % check if base case power flow converged
        if ~strcmp(contingencyCases{i,k}.mpcc.cpf.done_msg,'Base case power flow did not converge.')
            
            Nb = size(mpcb.bus,1);
            
            iter = mpcc.cpf.iterations;
            
            %% extract data
            V = mpcc.cpf.V;
            
            Pbase = sum(mpcb.bus(:,PD));
            Pdiff = sum(mpct.bus(:,PD) - mpcb.bus(:,PD));
            Pscale = Pbase + mpcc.cpf.lam * Pdiff;
            if CPFOptions.loadIterations % base load may be different, add it to total loadability
                pMax(i,k) = mpcc.cpf.max_lam * Pdiff + Pbase;
            else % base load is always the same, use this as zero point
                pMax(i,k) = mpcc.cpf.max_lam * Pdiff;
            end
            max_lam(i,k) = mpcc.cpf.max_lam;
            
            % power injections must be calculated from voltages and bus matrix
            Ybus = mpcc.cpf.Ybus;
            Sinj = zeros(size(mpcb.bus,1),iter+1);
            Pgen = zeros(size(mpcb.bus,1),iter+1);
            Qgen = zeros(size(Pgen));
            for ii=1:iter+1
                Sinj(:,ii) = V(:,ii).*conj(Ybus*V(:,ii))*mpcc.baseMVA;
                Pgen(:,ii) = mpcb.bus(:,PD) + mpcc.cpf.lam(ii)*(mpct.bus(:,PD)-mpcb.bus(:,PD)) ...
                    + real(Sinj(:,ii));
                Qgen(:,ii) = mpcb.bus(:,QD) + mpcc.cpf.lam(ii)*(mpct.bus(:,QD)-mpcb.bus(:,QD)) ...
                    + imag(Sinj(:,ii));
            end
            
            %% security constraints
            % check if voltage limits are satisfied
            ii = 0;
            s = 1;
            while s && ii <= iter;
                if size(V,2) < ii+1
                    disp('hej')
                end
                Vtmp = abs(V(:,ii+1));
                c1 = mpcb.bus(:,VMIN) >= Vtmp;
                c2 = mpcb.bus(:,VMAX) <= Vtmp;
                cond = sum(c1) + sum(c2);
                if cond == 0
                    ii = ii + 1;
                else
                    s = 0;
                    % save violation information
                    minLimitBuses = find(c1);
                    maxLimitBuses = find(c2);
                    if isempty(minLimitBuses) && ~isempty(maxLimitBuses)
                        securityLimitType(i,k) = 2; % max limit reached
                    elseif ~isempty(minLimitBuses) && isempty(maxLimitBuses)
                        securityLimitType(i,k) = 1; % min limit reached
                    else
                        error('Both max and min voltage limits violated'); % both max and min limits violated
                    end
                end
            end
            if ii > 0
                secure_lam(i,k) = mpcc.cpf.lam(ii);
                pSecure(i,k) = secure_lam(i,k) * Pdiff;
                secure(i,k) = 1;
            else
                secure_lam(i,k) = 0;
                pSecure(i,k) = 0;
            end
        else % base case power flow did not converge
            pMax(i,k) = 0;
            Sinj = [];
            Pgen = [];
            Qgen = [];
            Pscale = [];
            secure_lam(i,k) = 0;
            pSecure(i,k) = 0;
            secure(i,k) = 0;
            securityLimitType(i,k) = -1;
        end
        
        
        
        % store results in struct
        mpcc.cpf.secure_lam = secure_lam(i,k);
        mpcc.cpf.pMax = pMax(i,k);
        mpcc.cpf.pSecure = pSecure(i,k);
        mpcc.cpf.secure = secure(i,k);
        mpcc.cpf.Sinj = Sinj;
        mpcc.cpf.Pgen = Pgen;
        mpcc.cpf.Qgen = Qgen;
        mpcc.cpf.Pscale = Pscale;
        mpcc.cpf.securityLimitType = securityLimitType(i,k);
        mpcc.cpf.Pdiff = Pdiff;
        mpcc.cpf.Pbase = Pbase;
        
        contingencyCases{i,k}.mpcc = mpcc; % return struct with additional results
    end % loop over contingencies
end % loop of wind power

% collect description of contigencies in array
for i=1:CPFOptions.nContingencies
    s = strsplit(contingencyCases{i,1}.str,'\n');
    contingencies{i} = s{1};
end

res = struct();
res.pMax = pMax;
res.pSecure = pSecure;
res.max_lam = max_lam;
res.secure_lam = secure_lam;
res.secure = secure;
res.contingencies = contingencies;
res.pWind = CPFOptions.pWind;
res.securityLimitType = securityLimitType;

res.pMaxNminus1 = min(pMax);
res.pSecureNminus1 = min(pSecure);
res.Nminus1 = min(secure);



