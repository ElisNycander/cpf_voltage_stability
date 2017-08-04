% [res,contingencyCases] = processCPF(contingencyCases,CPFOptions)
% Collect and Plot results from cpfs contained in cell array contingencyCases
%
% The following extra results are stored in each mpcc.cpf:
% Sinj(bus x step) - net injected complex power at each bus
% Pgen(bus x step) - generated active power (injected power + load)
% Qgen(bus x step) - generated reactive power (injected power + load)
% Pmax             - system total load at maximum loadability point (0 if
%                    base case was not solvable)
% Psecure          - system total load at last point satisfying operating limits 
%

function [res,contingencyCases] = processCPF(contingencyCases,CPFOptions)

define_constants;

Ntot = length(contingencyCases);
Pmax = zeros(1,Ntot);
Psecure = zeros(1,Ntot);
max_lam = zeros(1,Ntot);
secure_lam = zeros(1,Ntot);
contingencies = cell(1,Ntot);
secure = zeros(1,Ntot);

% choose colormap
if ~isfield(CPFOptions.Plot,'colormap')
    ftmp = str2func('lines'); % default
else
    ftmp = str2func(CPFOptions.Plot.colormap);
end
line_colors = ftmp();
close;

for i=1:Ntot % loop over cases
    
    contingencies{i} = contingencyCases{i}.str;
    
    % check if base case power flow converged
    if ~strcmp(contingencyCases{i}.mpcc.cpf.done_msg,'Base case power flow did not converge.')
        
        
        mpcb = contingencyCases{i}.mpcb;
        mpct = contingencyCases{i}.mpct;
        mpcc = contingencyCases{i}.mpcc;
        
        Nb = size(mpcb.bus,1);
        
        iter = mpcc.cpf.iterations;
        
        %% extract data
        V = mpcc.cpf.V;
        
        Pbase = sum(mpcb.bus(:,PD));
        Pdiff = sum(mpct.bus(:,PD) - mpcb.bus(:,PD));
        Paxis = Pbase + mpcc.cpf.lam * Pdiff;
        Pmax(i) = mpcc.cpf.max_lam * Pdiff;
        max_lam(i) = mpcc.cpf.max_lam;
        
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
        while s
            Vtmp = abs(V(:,ii+1));
            cond = sum(mpcb.bus(:,VMIN) <= Vtmp) + ...
                sum(mpcb.bus(:,VMAX) >= Vtmp);
            if cond == 2*Nb
                ii = ii + 1;
            else
                s = 0;
            end
        end
        if ii > 0
            secure_lam(i) = mpcc.cpf.lam(ii);
            Psecure(i) = secure_lam(i) * Pdiff;
            secure(i) = 1;
        else
            secure_lam(i) = 0;
            Psecure(i) = 0;
        end
        
        
        
        %% Plot preparation
        if ~CPFOptions.Plot.disable
            % buses for which to Plot the voltage
            if isfield(CPFOptions.Plot,'voltageBuses')
                plot_v_bus = [];
                for ii = 1:length(CPFOptions.Plot.voltageBuses)
                    plot_v_bus = [plot_v_bus find( mpcb.bus(:,1) == CPFOptions.Plot.voltageBuses(ii))];
                end
            else
                plot_v_bus = 1:size(mpcb.bus,1);
            end
            Vplot = V(plot_v_bus,:);
            % create legend
            Nb_plot = length(plot_v_bus);
            vbus_legend = cell(1,Nb_plot);
            for ii=1:Nb_plot
                vbus_legend{ii} = ['Bus ' num2str(mpcb.bus(plot_v_bus(ii),1))];
            end
            vplot_colors = line_colors(plot_v_bus,:);
            
            % buses for which to Plot q generation
            if isfield(CPFOptions.Plot,'powerBuses')
                plot_q_bus = [];
                for ii=1:length(CPFOptions.Plot.powerBuses)
                    plot_q_bus = [plot_q_bus find( mpcb.bus(:,1) == CPFOptions.Plot.powerBuses(ii))];
                end
            else
                plot_q_bus = 1:size(mpcb.bus,1);
            end
            % create legend
            Nb_plot = length(plot_q_bus);
            qbus_legend = cell(1,Nb_plot);
            for ii=1:Nb_plot
                qbus_legend{ii} = ['Bus ' num2str(mpcb.bus(plot_q_bus(ii),1))];
            end
            qplot_colors = line_colors(plot_q_bus,:);
            
            
            if strcmp(CPFOptions.Plot.powerInjectedOrGenerated,'gen')
                Qplot = Qgen(plot_q_bus,:);
                Pplot = Pgen(plot_q_bus,:);
                pq_y_label = 'Q_{gen} (MW)';
                pp_y_label = 'P_{gen} (MW)';
            elseif strcmp(CPFOptions.Plot.powerInjectedOrGenerated,'inj')
                Qplot = imag(Sinj(plot_q_bus,:));
                Pplot = real(Sinj(plot_q_bus,:));
                pq_y_label = 'Q_{inj} (MW)';
                pp_y_label = 'P_{inj} (MW)';
            else
                error(['Unknown option for options.Plot.powerInjectedOrGenerated: ' CPFOptions.Plot.powerInjectedOrGenerated]);
            end
            
            %% plots
            if ~CPFOptions.Plot.subfigures
                if CPFOptions.Plot.pvCurve
                    % PV Plot
                    figure; pvh =  Plot(Paxis,abs(Vplot)); grid on;
                    legend(vbus_legend);
                    xlabel('P (MW)');
                    ylabel('V (p.u.)');
                    title(sprintf(['PV Plot\n' contingencyCases{i}.str]));
                end
                % PQ Plot
                if CPFOptions.Plot.pqCurve
                    figure; pqh = Plot(Paxis,Qplot); grid on;
                    legend(qbus_legend);
                    xlabel('P (MW)');
                    ylabel(pq_y_label);
                    title(sprintf(['PQ Plot\n' contingencyCases{i}.str]));
                end
                if CPFOptions.Plot.ppCurve
                    figure; pph = plot(Paxis,Pplot); grid on;
                    legend(qbus_legend);
                    xlabel('P (MW)');
                    ylabel(pp_y_label);
                    title(sprintf(['PP Plot\n' contingencyCases{i}.str]));
                end
            else % if subfigures
                figure;
                suptitle(contingencyCases{i}.str);
                subplot(2,2,1);
                
                if CPFOptions.Plot.pvCurve
                    % PV Plot
                    pvh =  plot(Paxis,abs(Vplot)); grid on;
                    legend(vbus_legend);
                    xlabel('P (MW)');
                    ylabel('V (p.u.)');
                    %title(sprintf(['PV Plot\n']));
                end
                
                %subplot(2,2,2);
                
                subplot(2,2,3);
                if CPFOptions.Plot.ppCurve % PP Plot
                    pph = plot(Paxis,Pplot); grid on;
                    xlabel('P (MW)');
                    ylabel(pp_y_label);
                    %title(sprintf(['PP Plot\n']));
                end
                
                subplot(2,2,4);
                if CPFOptions.Plot.pqCurve % PQ Plot
                    pqh = plot(Paxis,Qplot); grid on;
                    xlabel('P (MW)');
                    ylabel(pq_y_label);
                    %title(sprintf(['PQ Plot\n']));
                end
            end
            if CPFOptions.Plot.pvCurve
                for ii=1:length(pvh)
                    pvh(ii).Color = vplot_colors(ii,:);
                end
            end
            if CPFOptions.Plot.pqCurve
                for ii=1:length(pqh)
                    pqh(ii).Color = qplot_colors(ii,:);
                end
            end
            if CPFOptions.Plot.ppCurve
                for ii=1:length(pph)
                    pph(ii).Color = qplot_colors(ii,:);
                end
            end
        end
    else % base case power flow did not converge
        Pmax(i) = 0;
        Sinj = [];
        Pgen = [];
        Qgen = [];
        secure_lam(i) = 0;
        Psecure(i) = 0;
        secure(i) = 0;
    end
    
    % store results in struct
    mpcc.cpf.secure_lam = secure_lam(i);
    mpcc.cpf.Pmax = Pmax(i);
    mpcc.cpf.Psecure = Psecure(i);
    mpcc.cpf.secure = secure(i);
    mpcc.cpf.Sinj = Sinj;
    mpcc.cpf.Pgen = Pgen;
    mpcc.cpf.Qgen = Qgen;

    contingencyCases{i}.mpcc = mpcc; % return struct with additional results
end

res = struct();
res.Pmax = Pmax;
res.Psecure = Psecure;
res.max_lam = max_lam;
res.secure_lam = secure_lam;
res.secure = secure;
res.contingencies = contingencies;

res.pMaxNminus1 = min(Pmax);
res.pSecureNminus1 = min(Psecure);
res.Nminus1 = min(secure);


if ~CPFOptions.Plot.disable
    % Stability limits bar Plot
    if CPFOptions.Plot.voltageStabilityBar
        figure;
        Pdiff = Pmax - Psecure;
        bar(categorical(contingencies),[Psecure;Pdiff]','stacked');
        ylabel('P (MW)');
        set(gca,'xticklabel',contingencies);
        xticklabel_rotate([],30,[],'Fontsize',8)
        grid on;
        title('Voltage stability limits');
        % adjust y label
        ax = gca;
        ax.YLabel.Position = [-0.13 0.4 0];
    end
end



