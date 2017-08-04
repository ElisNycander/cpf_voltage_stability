% plot results in CPFResults from running CPFs with CPFOptions
%
%
%
%
function [] = plotCPF(CPFResults,CPFOptions)

shade_color = [0.5 1 0.8];
% define variables for shorter names
Pmax = CPFResults.pMax;
Pmax_min = CPFResults.pMaxNminus1;
Psecure = CPFResults.pSecure;
Psecure_min = CPFResults.pSecureNminus1;
secure_min = CPFResults.Nminus1;
Pdiff = Pmax - Psecure;
Pwind = CPFResults.pWind;
contingencies = CPFResults.contingencies;


if isfield(CPFResults,'pWind') % plot stability limits as function of pWind
    nPoints = length(CPFResults.pWind);
    if CPFOptions.Plot.separateContingencies % plot figure for each contingency
        for i=1:CPFOptions.nContingencies
            figure;
            if nPoints == 1
                bar([0 Pwind],[Psecure(i,:) Psecure(i,:);Pdiff(i,:) Pdiff(i,:)]','stacked');
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
else
end



