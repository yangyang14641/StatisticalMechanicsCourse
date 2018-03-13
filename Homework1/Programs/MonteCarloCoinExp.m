%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name:     MonteCarloCoinExp.m
%%%% Program purpose:  Monte Carlo Simulation of Coin Experiment
%%%% Program Aurthor:  Yang Yang
%%%% Creat Date:       25/2/2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters and Variables
N = 2^20;      % Numbers of experiments

%% Srochastic variables generation
% One coin situation
xi = zeros(1,N);                 % stochastic variables obey unifrom distribution
parfor i = 1:N
    xi(i) = rand();
end
eta = arrayfun(@Coin_Fun,xi);    % stochastic variables represent coin state
                                 % '1' for upside '0' for downside 
clear xi

% Two coin situation
xi = zeros(1,N);                 % stochastic variables obey unifrom distribution
parfor i = 1:N
    xi1(i) = rand();
    xi2(i) = rand();
end
eta1 = arrayfun(@Coin_Fun,xi1);    
eta2 = arrayfun(@Coin_Fun,xi2);

clear xi1 xi2

%% Statistics 
% only one coin
counter0 = 0;
counter1 = 0;
Freq_arry0 = zeros(1,N);
Freq_arry1 = zeros(1,N);
for i = 1:N
     if eta(i) == 0
         counter0 = counter0 + 1;
     else
         counter1 = counter1 + 1;
     end
     
     Freq_arry0(i) = counter0/i;
     Freq_arry1(i) = counter1/i;
end

% two coin situation
counter00 = 0;
counter01 = 0;
counter10 = 0;
counter11 = 0;

Freq_arry00 = zeros(1,N);
Freq_arry01 = zeros(1,N);
Freq_arry10 = zeros(1,N);
Freq_arry11 = zeros(1,N);
for i = 1:N
     if eta1(i) == 0 && eta2(i) == 0
         counter00 = counter00 + 1;
     elseif eta1(i) == 0 && eta2(i) == 1
         counter01 = counter01 + 1;
     elseif eta1(i) == 1 && eta2(i) == 0
         counter10 = counter10 + 1;    
     else
         counter11 = counter11 + 1;
     end
     
     Freq_arry00(i) = counter00/i;
     Freq_arry01(i) = counter01/i;
     Freq_arry10(i) = counter10/i;
     Freq_arry11(i) = counter11/i;     
end



%% plot results
% one coin
fig1 = figure(1);
set(fig1,'Color',[1 1 1]);
subplot(1,2,1);
semilogx(1/2*ones(1,N*10),'k--','linewidth',2);
grid on; hold on;
semilogx(Freq_arry1,'r','linewidth',2);
legend('Probability','Frequency');
title('Upside');
xlabel('Expriment times');
ylabel('Frequency & Probability');
hold off;

subplot(1,2,2)
semilogx(1/2*ones(1,N*10),'k--','linewidth',2);
grid on; hold on;
semilogx(Freq_arry0,'b','linewidth',2);
legend('Probability','Frequency');
title('Downside');
xlabel('Expriment times');
ylabel('Frequency & Probability');
hold off;

% two coin
fig2 = figure(2);
set(fig2,'Color',[1 1 1]);
subplot(2,2,1)
semilogx(1/4*ones(1,N*10),'k--','linewidth',2);
grid on; hold on;
semilogx(Freq_arry11,'r','linewidth',2);
legend('Probability','Frequency');
title('Upside Upside');
xlabel('Expriment times');
ylabel('Frequency & Probability');
hold off;

subplot(2,2,2)
semilogx(1/4*ones(1,N*10),'k--','linewidth',2);
grid on; hold on;
semilogx(Freq_arry10,'b','linewidth',2);
legend('Probability','Frequency');
title('Upside Downside');
xlabel('Expriment times');
ylabel('Frequency & Probability');
hold off;

subplot(2,2,3)
semilogx(1/4*ones(1,N*10),'k--','linewidth',2);
grid on; hold on;
semilogx(Freq_arry01,'g','linewidth',2);
legend('Probability','Frequency');
title('Downside Upside');
xlabel('Expriment times');
ylabel('Frequency & Probability');
hold off;

subplot(2,2,4)
semilogx(1/4*ones(1,N*10),'k--','linewidth',2);
grid on; hold on;
semilogx(Freq_arry00,'m','linewidth',2);
legend('Probability','Frequency');
title('Downside Downside');
xlabel('Expriment times');
ylabel('Frequency & Probability');
hold off;

clear