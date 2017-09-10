% plot results in CPFProcessedResults from running CPFs with CPFOptions
%
%
%
%
function [] = plotCPF(CPFProcessedResults,CPFResults,CPFOptions)


% choose colormap
if ~isfield(CPFOptions.Plot,'colormap')
    plotLineColors = lineColors('lines');
else
    plotLineColors = lineColors(CPFOptions.Plot.colormap);
end

shade_color = [0.5 1 0.8];
% define variables for shorter names
Pmax = CPFProcessedResults.pMax;
Pmax_min = CPFProcessedResults.pMaxNminus1;
Psecure = CPFProcessedResults.pSecure;
Psecure_min = CPFProcessedResults.pSecureNminus1;
secure_min = CPFProcessedResults.Nminus1;
Pdiff = Pmax - Psecure;
Pwind = CPFProcessedResults.pWind;
contingencies = CPFProcessedResults.contingencies;


if CPFOptions.Plot.stabilityMargin % plot stability limits as function of pWind
    nPoints = length(CPFProcessedResults.pWind);
    if CPFOptions.nWindPoints == 1 % plot all contingencies in same figure
        
        figure;
        if CPFOptions.nContingencies == 1
            tempCell = {'N/A',contingencies{1}};
            bar(categorical(tempCell),[0 Psecure(:);0 Pdiff(:)]','stacked');
        else
            tempCell = contingencies;
            bar(categorical(contingencies),[Psecure(:) Pdiff(:)],'stacked');
        end
        ylabel('P (MW)');
        set(gca,'xticklabel',tempCell);
        try % xticklabel sometimes gives error, not certain why 
            xticklabel_rotate([],30,[],'Fontsize',8)
        catch
        end
        grid on;
        title(sprintf('Voltage stability limits, %i MW wind',Pwind));
        % adjust y label
        ax = gca;
        ax.YLabel.Position = [-0.13 0.4 0];
        
    elseif CPFOptions.Plot.separateContingencies % plot figure for each contingency
        for i=1:CPFOptions.nContingencies
            figure;
            if nPoints == 1
                bar([0 Pwind],[0 Psecure(i,:);0 Pdiff(i,:)]','stacked');
            else
                bar(Pwind,[Psecure(i,:);Pdiff(i,:)]','stacked');
            end
            ylabel('P (MW)');
            xlabel('P_{wind} (MW)');
            grid on;
            title(sprintf(['Voltage stability margin\n' contingencies{i}]));
        end
    else % plot curves for all contingencies together
        
        % maximum loadability limits
        figure;
        ph1 = plot(Pwind,Pmax,'-o');
        grid on; hold on;
        xlabel('P_{wind} (MW)');
        ylabel('P_{max} (MW)');
        title('Loadability margin');
        legend(contingencies);
        % shade area
        a1 = area(Pwind,Pmax_min);
        a1.BaseLine.LineStyle = 'none';
        a1.FaceColor = shade_color;
        a1.LineStyle = 'none';
        alpha(a1,0.2);
        
        % security limits
        figure;
        ph2 = plot(Pwind,Psecure,'-o');
        grid on; hold on;
        xlabel('P_{wind} (MW)');
        ylabel('P_{secure} (MW)');
        title('Security margin');
        legend(contingencies);
        % shade area
        a2 = area(Pwind,Psecure_min);
        a2.BaseLine.LineStyle = 'none';
        a2.FaceColor = shade_color;
        a2.LineStyle = 'none';
        alpha(a2,0.2);
        
        % add security index
        bh = bar(Pwind,50*secure_min);
        alpha(bh,0.9);
        bh.BarWidth = 0.2;
        bh.LineStyle = '-';
        bh.FaceColor = [1 1 1];
    end
end

if CPFOptions.Plot.individualCases
    
    for k=1:size(CPFResults,1)
        for l=1:size(CPFResults,2)
            
            if ~strcmp(CPFResults{k,l}.mpcc.cpf.done_msg,'Base case power flow did not converge.')
                mpcb = CPFResults{k,1}.mpcb;
                cpf = CPFResults{k,l}.mpcc.cpf;
                caseString = CPFResults{k,l}.str;
                
                %% Plot preparation
                % buses for which to Plot the voltage
                if isfield(CPFOptions.Plot,'voltageBuses')
                    plot_v_bus = [];
                    for ii = 1:length(CPFOptions.Plot.voltageBuses)
                        plot_v_bus = [plot_v_bus find( mpcb.bus(:,1) == CPFOptions.Plot.voltageBuses(ii))];
                    end
                end
                if isempty(plot_v_bus)
                    plot_v_bus = 1:size(mpcb.bus,1);
                end
                Vplot = cpf.V(plot_v_bus,:);
                % create legend
                Nb_plot = length(plot_v_bus);
                vbus_legend = cell(1,Nb_plot);
                for ii=1:Nb_plot
                    vbus_legend{ii} = ['Bus ' num2str(mpcb.bus(plot_v_bus(ii),1))];
                end
                vplot_colors = plotLineColors(plot_v_bus,:);
                
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
                qplot_colors = plotLineColors(plot_q_bus,:);
                
                
                if strcmp(CPFOptions.Plot.powerInjectedOrGenerated,'gen')
                    Qplot = cpf.Qgen(plot_q_bus,:);
                    Pplot = cpf.Pgen(plot_q_bus,:);
                    pq_y_label = 'Q_{gen} (MW)';
                    pp_y_label = 'P_{gen} (MW)';
                elseif strcmp(CPFOptions.Plot.powerInjectedOrGenerated,'inj')
                    Qplot = imag(cpf.Sinj(plot_q_bus,:));
                    Pplot = real(cpf.Sinj(plot_q_bus,:));
                    pq_y_label = 'Q_{inj} (MW)';
                    pp_y_label = 'P_{inj} (MW)';
                else
                    error(['Unknown option for options.Plot.powerInjectedOrGenerated: ' CPFOptions.Plot.powerInjectedOrGenerated]);
                end
                
                %% plots
                if ~CPFOptions.Plot.subfigures
                    if CPFOptions.Plot.pvCurve
                        % PV Plot
                        figure; pvh =  plot(cpf.Pscale,abs(Vplot)); grid on;
                        legend(vbus_legend);
                        xlabel('P (MW)');
                        ylabel('V (p.u.)');
                        title(sprintf(['Bus voltages\n' CPFResults{k,l}.str]));
                    end
                    % PQ Plot
                    if CPFOptions.Plot.pqCurve
                        figure; pqh = plot(cpf.Pscale,Qplot); grid on;
                        legend(qbus_legend);
                        xlabel('P (MW)');
                        ylabel(pq_y_label);
                        title(sprintf(['Reactive power\n' CPFResults{k,l}.str]));
                    end
                    if CPFOptions.Plot.ppCurve
                        figure; pph = plot(cpf.Pscale,Pplot); grid on;
                        legend(qbus_legend);
                        xlabel('P (MW)');
                        ylabel(pp_y_label);
                        title(sprintf(['Active power\n' CPFResults{k,l}.str]));
                    end
                else % if subfigures
                    figure;
                    suptitle(CPFResults{k,l}.str);
                    subplot(2,2,1);
                    
                    if CPFOptions.Plot.pvCurve
                        % PV Plot
                        pvh =  plot(cpf.Pscale,abs(Vplot)); grid on;
                        legend(vbus_legend);
                        xlabel('P (MW)');
                        ylabel('V (p.u.)');
                        %title(sprintf(['PV Plot\n']));
                    end
                    
                    %subplot(2,2,2);
                    
                    subplot(2,2,3);
                    if CPFOptions.Plot.ppCurve % PP Plot
                        pph = plot(cpf.Pscale,Pplot); grid on;
                        xlabel('P (MW)');
                        ylabel(pp_y_label);
                        %title(sprintf(['PP Plot\n']));
                    end
                    
                    subplot(2,2,4);
                    if CPFOptions.Plot.pqCurve % PQ Plot
                        pqh = plot(cpf.Pscale,Qplot); grid on;
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
            end % if CPF converged
        end % for columns
    end % for rows
end % if plot individual plots



