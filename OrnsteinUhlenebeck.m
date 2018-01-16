%%%% Author: Belén Sancristóbal
%%%% Date: 1st April 2016
%%%% Description: Analysis of the PSD of an Ornstein-Uhlenbeck process
%%%% according to different "baseline" activities.
%%%% Comment on "Characterization of Subthreshold Voltage Fluctuations in
%%%% Neuronal Membranes", by M. Rudolph and A. Destexhe
%%%% B. Lindner and A.Longtin, Neural Computation, 2006

clearvars 

%% Ornstein-Uhlenebeck process (eq.2.3/2.4)
%% dOU/dt = -(OU-OU_o)/tau + sqrt(2*sd^2/tau)*epsilon(t)
%% epsilon(t) is an independent gaussian white noise

tau = 10.^linspace(log10(0.01),log10(10),50); %s
var = 1;
OU_o = 0;
dt = 0.001*ones(1,length(tau));%tau*10e-4; %s

R = 1; %number of realizations
NTrials = 30; %number of trials/pseudo-trials
TF = 5; %end epoch
T0 = -1.5; %start epoch

time_corr = 0.2; %s

% corr_OU = zeros(length(tau),R*NTrials);
% corr_REST = zeros(length(tau),R*NTrials);

corr_OU_high = zeros(1,length(tau));
corr_OU_low = zeros(1,length(tau));

numeric_var = zeros(1,length(tau));

AC_OU = cell(1,length(tau));
AC_OU_ZeroLag = zeros(1,length(tau));
AC_fft_OU = cell(1,length(tau));

tic;

for m = 1:length(tau)
    tau(m)
    
    t = (0:dt(m):10000);
    t_trial = (T0:dt(m):TF);

    LI = t(end)-TF-abs(T0);  %length of interval
    d = TF + abs(T0);   %minimum distance
    E = LI - (NTrials-1)*d;  %excess space for points
    %generate N+1 random values; 
    Ro = rand(NTrials+1,1);     %random vector
    %normalize so that the extra space is consumed
    % extra value is the amount of extra space "unused"
    Rn = E*Ro(1:NTrials)/sum(Ro); %normalize
    %spacing of points
    S=d*ones(NTrials,1)+Rn; 

    %location of points, adjusted to "start" at 0
    P = cumsum(S)-1;
    P(P > t(end)-TF) = [];
    % PSEUDO-TRIALS
    Trigger = round(P/dt(m))+1;
    % TRIALS
    locEF = Trigger + round(TF/dt(m));

    % BASELINE
    locB0 = Trigger + round(T0/dt(m));

    TS = 3; %s
    loom = intensity((0:dt(m):TS));
    lmax = 0.005;
    lmin = 0.001;
    loom = ((lmax-lmin)/(max(loom)-min(loom)))*(loom-min(loom))+lmin;

    [~,tloc_ini] = min(abs(t_trial-0));
    [~,tloc_end] = min(abs(t_trial-TS));
    
    OU = zeros(R,length(t));  %%Ornstein-Uhlenebeck
    REST = zeros(R,length(t));

    REST(:,1) = randn(R,1);
    OU(:,1) = randn(R,1);
    
    for i = 2:length(t)-1
        if ~isempty(find(i >= Trigger & i <= Trigger+round(TS/dt(m)), 1))
           I = loom(i-Trigger(i >= Trigger & i <= Trigger+round(TS/dt(m)))+1); 
        else 
           I = 0;
        end
        REST(:,i) = REST(:,i-1) - dt(m)*(REST(:,i-1)-repmat(OU_o,R,1))/tau(m) + sqrt(2*dt(m)*var/tau(m))*randn(R,1);
        OU(:,i) = OU(:,i-1) - dt(m)*(OU(:,i-1)-repmat(OU_o,R,1))/tau(m) + sqrt(2*dt(m)*var/tau(m))*randn(R,1) + I*ones(R,1);
    end
 
    numeric_var(m) = mean(std(REST,[],2).^2,1);
    
    OU_zscore = (OU - repmat(mean(OU,2),1,size(t,2)))./repmat(std(OU,[],2),1,size(t,2)); %% If Statistics Toolbox zscore(psd_bands,[],2);
    REST_zscore = (REST - repmat(mean(REST,2),1,size(t,2)))./repmat(std(REST,[],2),1,size(t,2)); %% If Statistics Toolbox zscore(psd_bands,[],2);
    
    trials_OU = [];
    pseudo_REST = [];

    for k = 1:R
        trials_OU = cat(2,trials_OU,arrayfun(@(i) squeeze(OU_zscore(k,locB0(i):locEF(i)))', (1:length(locB0)), 'UniformOutput', false));
        pseudo_REST = cat(2,pseudo_REST,arrayfun(@(i) squeeze(REST_zscore(k,locB0(i):locEF(i)))', (1:length(locB0)), 'UniformOutput', false));
    end

    figure
    subplot(2,3,[1 3])
    plot(t_trial,mean(cell2mat(trials_OU),2))
    hold all
    plot(t_trial,mean(cell2mat(pseudo_REST),2),'k')
    xlim([T0 TF])
    xlabel('time (s)')
    legend({'OU','REST'})
    title(['tau = ',num2str(round(tau(m)*10^2)/10^2),'s'])
    set(findall(gcf,'-property','FontSize'),'FontSize',20)

    %%% OU %%%
    bsln_trials = cell2mat(cellfun(@(w) mean(w(1:tloc_ini-1)), trials_OU, 'UniformOutput', false));
    median_bsln = median(bsln_trials);

    index_high = find(bsln_trials >= median_bsln);
    index_low = find(bsln_trials < median_bsln);

    mean_high_OU = mean(cell2mat(trials_OU(index_high)),2)';
    sd_high_OU = std(cell2mat(trials_OU(index_high)),[],2)';
    mean_low_OU = mean(cell2mat(trials_OU(index_low)),2)';
    sd_low_OU = std(cell2mat(trials_OU(index_low)),[],2)';

    subplot(2,3,4)
    plot(t_trial,mean_high_OU)
    hold all
    plot(t_trial,mean_low_OU)
    xlim([T0 TF])
    xlabel('time (s)')
    ylabel('OU')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)


    %%% REST %%%
    bsln_trials = cell2mat(cellfun(@(w) mean(w(1:tloc_ini-1)), pseudo_REST, 'UniformOutput', false));
    median_bsln = median(bsln_trials);

    index_high = find(bsln_trials >= median_bsln);
    index_low = find(bsln_trials < median_bsln);

    mean_high_REST = mean(cell2mat(pseudo_REST(index_high)),2)';
    sd_high_REST = std(cell2mat(pseudo_REST(index_high)),[],2)';
    mean_low_REST = mean(cell2mat(pseudo_REST(index_low)),2)';
    sd_low_REST = std(cell2mat(pseudo_REST(index_low)),[],2)';

    subplot(2,3,5)
    plot(t_trial,mean_high_REST)
    hold all
    plot(t_trial,mean_low_REST)
    xlim([T0 TF])
    xlabel('time (s)')
    ylabel('REST')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)

    subplot(2,3,6)
    plot(t_trial,mean_high_OU-mean_high_REST)
    hold all
    plot(t_trial,mean_low_OU-mean_low_REST)
    xlim([T0 TF])
    xlabel('time (s)')
    ylabel('OU-REST')
    set(findall(gcf,'-property','FontSize'),'FontSize',20)
    
    %% CORRELATIONS WITH SOUND %%
    %% trial-averaged based
    corr_OU_high(m) = corr((mean_high_OU(tloc_ini:tloc_end)-mean_high_REST(tloc_ini:tloc_end))',loom','type','Spearman');
    corr_OU_low(m) = corr((mean_low_OU(tloc_ini:tloc_end)-mean_low_REST(tloc_ini:tloc_end))',loom','type','Spearman');    
    
    AC_OU{m} = mean(cell2mat(arrayfun(@(m) xcov(detrend(REST(m,:),'constant'),'unbiased'), (1:R), 'UniformOutput', false)'),1);
    AC_OU_ZeroLag(m) = AC_OU{m}(length(t));
end


figure
plot(tau,corr_OU_high,'-o')
hold all
plot(tau,corr_OU_low,'-o')
xlabel('\tau (s)')
ylabel('CC')
legend({'High';'Low'})
set(findall(gcf,'-property','FontSize'),'FontSize',20)

figure
plot(numeric_var,AC_OU_ZeroLag,'o')
hold all
plot(linspace(0.7,1.5,10),linspace(0.7,1.5,10),'--','color','k')
xlabel('\sigma^2')
ylabel('ACF (0)')
set(findall(gcf,'-property','FontSize'),'FontSize',20)


toc;