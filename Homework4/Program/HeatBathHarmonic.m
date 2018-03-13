%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name : HeatBathHarmonic.m
%%%% Program Author: Yang Yang
%%%% Data : 12/3/2016
%%%% Version: V1.0
%%%% copyright: all  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation parameters
global Ms a b
% Ms = 0.1;
% Ms = 1;
Ms = 5;         
a = 1;
b = 1;

y0 = [1 0 0];
tspan = [0 250];
%% Computing
options = odeset('RelTol',1e-12);
[t,y] = ode45(@HBHarmonic,tspan,y0,options);
p = y(:,1);
q = y(:,2);
eta = y(:,3);


[t_hat,y_hat] = ode45(@Harmonic,tspan,y0(1:2),options);
p_hat = y_hat(:,1);
q_hat = y_hat(:,2);

%% Plot
f1 = figure(1);
set(f1,'Color',[1 1 1 ]);
set(f1,'Name','Harmonic Oscillator phase plane');
plot(p,q,'r',p_hat,q_hat,'b','linewidth',2);
xlabel('p');
ylabel('q');
title('Oscillators');
legend('Heat bath Oscillator','Harmonic Oscillator');

f2 = figure(2);
set(f2,'Color',[1 1 1 ]);
set(f2,'Name','Harmonic Oscillator time series');
subplot(4,1,1)
plot(t,p);
xlabel('t');
ylabel('p');
title('Heat bath Oscillator p(t)');
subplot(4,1,2)
plot(t,q);
xlabel('t');
ylabel('q');
title('Heat bath Oscillator q(t)');
subplot(4,1,3)
plot(t_hat,p_hat);
xlabel('t');
ylabel('p');
title('Harmonic Oscillator Oscillator p(t)');
subplot(4,1,4)
plot(t_hat,q_hat);
xlabel('t');
ylabel('q');
title('Harmonic Oscillator Oscillator q(t)');




%% Program End