%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name : MDDataProcess.m
%%%% Program Author: Yang Yang
%%%% Data : 12/3/2016
%%%% Version: V1.0
%%%% copyright: all
%%%% input file format: result.txt    
%%%% data column (t, x1, y1, u1, v1, k1, p1, x2, y2, u2, v2, k2, p2, x3, y3, u3, v3, k3, p3)
%%%% Thanks to our programmer: Huang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data input module
filename = 'result.txt';        % data file name goes here
Data = importdata(filename);    % import command
[M,N] = size(Data);             % find dimension of data


%% Regroup data
% time array
t = Data(:,1);

% x-coordiante arrays
NX = (N-1)/6;
X = zeros(M,NX);
for j = 1:NX
   X(:,j) = Data(:,2 + (j-1)*6);
end

% y-coordinate arrays
NY = (N-1)/6;
Y = zeros(M,NX);
for j = 1:NY
    Y(:,j) = Data(:,3 + (j-1)*6);
end

% Velocity x component arrays
NVX = (N-1)/6;
VX = zeros(M,NVX);
for j = 1:NVX
   VX(:,j) = Data(:,4 + (j-1)*6);
end

% Velocity y component arrays
NVY = (N-1)/6;
VY = zeros(M,NVY);
for j = 1:NVY
    VY(:,j) = Data(:,5 + (j-1)*6);
end

% knetic energy arrays
NKinetic = (N-1)/6;
Kinetic = zeros(M,NKinetic);
for j = 1:NKinetic
   Kinetic(:,j) = Data(:,6 + (j-1)*6);
end

% Potential energy arrays
NPotential = (N-1)/6;
Potential = zeros(M,NPotential);
for j = 1:NPotential
    Potential(:,j) = Data(:,7 + (j-1)*6);
end

% Particle's total energy arrays
NEnergy = (N-1)/6;
Energy = zeros(M,NEnergy);
for j = 1:NEnergy
    Energy(:,j) = Kinetic(:,j) + Potential(:,j);
end


% Ensemble total energy
NEnsembleEnergy = 1;
EnsembleEnergy = zeros(M,NEnsembleEnergy);
EnsembleEnergy = sum(Energy,2);

%% GUI plot and output module
% Particle Position trace Movie    $q_{i}$ 
f1 = figure(1);
set(f1,'Name','Trace movie');
set(f1,'Color',[1 1 1]);
subplot(3,1,1);
comet(X(:,1),Y(:,1));
box on
title('No.1 Particle Trace');
xlabel('x');
ylabel('y');
subplot(3,1,2);
comet(X(:,2),Y(:,2));
box on
title('No.2 Particle Trace');
xlabel('x');
ylabel('y');
subplot(3,1,3);
comet(X(:,3),Y(:,3));
box on
title('No.3 Particle Trace');
xlabel('x');
ylabel('y');

% Particle Trace   $\q_{i}$
f2 = figure(2);
set(f2,'Name','Particle Trace');
set(f2,'Color',[1 1 1]);
plot(X(:,1),Y(:,1),'r','linewidth',2);
hold on
plot(X(:,2),Y(:,2),'b','linewidth',2);
hold on
plot(X(:,3),Y(:,3),'m','linewidth',2);
hold on
legend('Particle No.1','Particle No.2','Particle No.3')
title('Particle Trace');
xlabel('x');
ylabel('y');
hold off


% Plot Energy  
f3 = figure(3);
set(f3,'Name','Energy Plot');
set(f3,'Color',[1 1 1]);
subplot(2,2,1);
plot(t,Energy(:,1));
title('Particle No.1 Energy');
xlabel('Time');
ylabel('Energy');
subplot(2,2,2);
plot(t,Energy(:,2));
title('Particle No.2 Energy');
xlabel('Time');
ylabel('Energy');
subplot(2,2,3);
plot(t,Energy(:,3));
title('Particle No.3 Energy');
xlabel('Time');
ylabel('Energy');
subplot(2,2,4);
plot(t,EnsembleEnergy);
title('Ensemble Energy');
xlabel('Time');
ylabel('Energy');



% Plot Phase space (Here we decompose the $3\times 4D = 12D$ Phase space into six 2D space)
f4 = figure(4);
set(f4,'Name','Phase space');
set(f4,'Color',[1 1 1]);
% $(\p_{1,x},q_{1,x})$ 
subplot(3,2,1);
plot(X(:,1),VX(:,1),'r')
title('Particle No.1:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{1,y},q_{1,y})$ 
subplot(3,2,2);
plot(Y(:,1),VY(:,1),'r')
title('Particle No.1:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{2,x},q_{2,x})$ 
subplot(3,2,3);
plot(X(:,2),VX(:,2),'b')
title('Particle No.2:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{2,y},q_{2,y})$ 
subplot(3,2,4);
plot(Y(:,2),VY(:,2),'b')
title('Particle No.2:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{3,x},q_{3,x})$ 
subplot(3,2,5);
plot(X(:,3),VX(:,3),'m')
title('Particle No.3:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{3,y},q_{3,y})$ 
subplot(3,2,6);
plot(Y(:,3),VY(:,3),'m')
title('Particle No.3:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
%% Program End