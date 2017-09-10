function plotPFSurface(filename,CPFOptions)

%xAxisUnit = 'power_angle'; % power_factor / power_angle

if nargin < 1
    load('pf_results.mat') % default file
else
    load(filename);
end

if nargin < 2
    CPFOptions.xScale = 'qp_ratio'; % Q/P ratio default scale
end

titlestring = strrep(filename,'_',' ');
titlestring = strrep(titlestring,'.mat','');

if strcmp(CPFOptions.xScale,'power_angle')
    qScale = 180/pi*results.powerAngle;
    qScaleLabel = sprintf('Power angle (%c)',char(176));
else
    qScale = tan(results.powerAngle);
    qScaleLabel = 'Q/P ratio';
end

pWind = results.pWind;
powerAngle = results.powerAngle;
loadabilityMargin = results.loadabilityMargin;
securityMargin = results.securityMargin;
securityLimitType = results.securityLimitType;

results.nosePLossFraction(results.nosePLossFraction == 0) = NaN;
results.basePLossFraction(results.basePLossFraction == 0) = NaN;
loadabilityMargin(loadabilityMargin == 0) = NaN;
securityMargin(securityMargin == 0) = NaN;


figure;
surf(results.pWind,qScale,loadabilityMargin);
colorbar;
title(sprintf(['Maximum loadability\n' titlestring]))
xlabel('Wind power (MW)');
ylabel(qScaleLabel);
zlabel('Maximum load (MW)');
figure;
surf(results.pWind,qScale,securityMargin);
colorbar;
title(sprintf(['Maximum load under voltage security\n' titlestring]))
xlabel('Wind power (MW)');
ylabel(qScaleLabel);
zlabel('Security margin (MW)');
figure;
surf(results.pWind,qScale,securityLimitType);
colorbar;
title(sprintf(['Security violation type\n' titlestring]))
xlabel('Wind power (MW)');
ylabel(qScaleLabel);
zlabel('Security violation type');
figure;
surf(results.pWind,qScale,results.basePLossFraction*100);
colorbar;
title(sprintf(['Base case active power losses\n' titlestring]))
xlabel('Wind power (MW)');
ylabel(qScaleLabel);
zlabel('Losses (%)');
figure;
surf(results.pWind,qScale,results.nosePLossFraction*100);
colorbar;
title(sprintf(['Nose point case active power losses\n' titlestring]))
xlabel('Wind power (MW)');
ylabel(qScaleLabel);
zlabel('Losses (%)');

% plot separate curves for different power factors
figure;
plot(results.pWind,results.loadabilityMargin');

npf = size(results.loadabilityMargin,1);
labels = cell(1,npf);
for i=1:npf
    if results.powerAngle(i) >= 0 
        str = 'lag';
    else
        str = 'lead';
    end
    %labels{i} = sprintf(['PF %0.3f ' str],cos(results.powerAngle(i))); 
    labels{i} = sprintf('Q/P = %0.2f',tan(results.powerAngle(i)));
end
legend(labels);
grid on;
title(sprintf(['Maximum loadability\n' titlestring]));
xlabel('Wind power (MW)');
ylabel('Total load (MW)');


%% optimal loadability

optimalLoadability = zeros(1,npf);
optimalWindPower = zeros(1,npf);
for i=1:npf
    [optimalLoadability(i),idx] = max(results.loadabilityMargin(i,:));
    optimalWindPower(i) = results.pWind(idx);
end

figure; % power angle scale
subplot(1,2,1);
plot(qScale,optimalWindPower);
xlabel(qScaleLabel);
grid on;
ylabel('Wind power (MW)');

subplot(1,2,2);
plot(qScale,optimalLoadability);
xlabel(qScaleLabel);
grid on;
ylabel('Maximum loadability (MW)');

suptitle(sprintf(['Optimal wind power for maximum loadability\n' titlestring]));

%% optimal loadability under voltage security

optimalSecureLoadability = zeros(1,npf);
optimalSecureWindPower = zeros(1,npf);
for i=1:npf
    [optimalSecureLoadability(i),idx] = max(results.securityMargin(i,:));
    optimalSecureWindPower(i) = results.pWind(idx);
end

figure; % power angle scale
subplot(1,2,1);
plot(qScale,optimalSecureWindPower);
xlabel(qScaleLabel);
grid on;
ylabel('Wind power (MW)');

subplot(1,2,2);
plot(qScale,optimalSecureLoadability);
xlabel(qScaleLabel);
grid on;
ylabel('Maximum load (MW)');

suptitle(sprintf(['Optimal wind power for maximum loadability under voltage security\n' titlestring]));

% figure; % power factor scale
% subplot(1,2,1);
% plot(results.powerFactor,optimalWindPower);
% xlabel('Power factor')
% 
% subplot(1,2,2);
% plot(results.powerFactor,optimalLoadability);
% xlabel('Power factor')
% subplot(1,2,1)
% grid on;
% ylabel('Optimal amount of wind power (MW)');
% subplot(1,2,2);
% grid on;
% ylabel('Maximum loadability margin (MW)');
% suptitle(titlestring);

%% optimal losses
figure;
plot(results.pWind,results.basePLossFraction');
grid on;
xlabel('Wind power (MW)')
ylabel('Losses (%)');
legend(labels)
title(sprintf(['Losses in base case (500 MW load)\n' titlestring]));

optimalLosses = zeros(1,npf);
optimalLossWindPower = zeros(1,npf);
optimalNoseLosses = zeros(1,npf);
optimalNoseLossWindPower = zeros(1,npf);

for i=1:npf
    [optimalLosses(i),idx] = min(results.basePLossFraction(i,:));
    optimalLossWindPower(i) = results.pWind(idx);
    
    [optimalNoseLosses(i),idx] = min(results.nosePLossFraction(i,:));
    optimalNoseLossWindPower(i) = results.pWind(idx);
end

figure;
subplot(1,2,1);
plot(qScale,optimalLossWindPower,qScale,optimalNoseLossWindPower);
xlabel(qScaleLabel);
ylabel('Wind power (MW)');
legend('Losses at base load','Losses at nose point');
grid on;
subplot(1,2,2);
plot(qScale,100*optimalLosses,qScale,100*optimalNoseLosses);
xlabel(qScaleLabel);
ylabel('Losses at base load (%)')
grid on;
suptitle(sprintf(['Optimal wind power for minimal losses\n' titlestring]));

%% optimal losses at nose point
% optimalNoseLosses = zeros(1,npf);
% optimalNoseLossWindPower = zeros(1,npf);
% 
% for i=1:npf
%     [optimalNoseLosses(i),idx] = min(results.nosePLossFraction(i,:));
%     optimalNoseLossWindPower(i) = results.pWind(idx);
% end
% 
% figure;
% subplot(1,2,1);
% plot(qScale,optimalNoseLossWindPower);
% xlabel(qScaleLabel);
% ylabel('Wind power (MW)');
% grid on;
% subplot(1,2,2);
% plot(qScale,100*optimalNoseLosses);
% xlabel(qScaleLabel);
% ylabel('Losses at CPF nose point (%)')
% grid on;
% suptitle(sprintf(['Optimal wind power for minimal losses at nose point\n' titlestring]));

%% common plot with all optima

figure;
plot(qScale,[optimalWindPower; optimalSecureWindPower; optimalLossWindPower; optimalNoseLossWindPower]);
xlabel(qScaleLabel);
ylabel('Wind power (MW)');
grid on;
title('Optimal wind power for different criteria');
legend('Max loadability','Max loadability under voltage security','Minimal losses for base load',...
    'Minimal losses for nose point');