%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name : MDDataAnalysis.m
%%%% Program Author: Yang Yang
%%%% Data : 30/4/2016
%%%% Version: V1.0
%%%% copyright: Author only!
%%%% input file format: result.txt    
%%%% data column ( 
%%%% t, x1, y1, u1, v1, fx1, fy1, k1, p1, 
%%%%    x2, y2, u2, v2, fx2, fy2, k2, p2, 
%%%% ......
%%%%    x10, y10, u10, v10, k10, p10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data input module
filename = 'result.txt';        % data file name goes here
Data = importdata(filename);    % import command
[M,N] = size(Data);             % find dimension of data

problemDim = 2;                 % problem dimension

%% Regroup data
col_num = 8;
index = 1;

% time array
t = Data(:,1);
index = index + 1;

% x-coordiante arrays
NX = (N-1)/col_num;
X = zeros(M,NX);
for j = 1:NX
   X(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;

% y-coordinate arrays
NY = (N-1)/col_num;
Y = zeros(M,NX);
for j = 1:NY
    Y(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;

% Velocity x component arrays
NVX = (N-1)/col_num;
VX = zeros(M,NVX);
for j = 1:NVX
   VX(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;

% Velocity y component arrays
NVY = (N-1)/col_num;
VY = zeros(M,NVY);
for j = 1:NVY
    VY(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;

% internal force x component arrays
NFX = (N-1)/col_num;
FX = zeros(M,NVY);
for j = 1:NVY
    FX(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;

% internal force y component arrays
NFY = (N-1)/col_num;
FY = zeros(M,NVY);
for j = 1:NVY
    FY(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;

% knetic energy arrays
NKinetic = (N-1)/col_num;
Kinetic = zeros(M,NKinetic);
for j = 1:NKinetic
   Kinetic(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;

% Potential energy arrays
NPotential = (N-1)/col_num;
Potential = zeros(M,NPotential);
for j = 1:NPotential
    Potential(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;


% Particle's total energy arrays
NEnergy = (N-1)/col_num;
Energy = zeros(M,NEnergy);
for j = 1:NEnergy
    Energy(:,j) = Kinetic(:,j) + Potential(:,j);
end


% Ensemble total energy

Ensemble_Kinetic = sum(Kinetic,2);
Ensemble_Potential = sum(Potential,2);
Ensemble_Energy = Ensemble_Kinetic + Ensemble_Potential;

%% plot data to get equlibrium time step
figure('Color',[1 1 1]);
plot(Ensemble_Kinetic);
title('Ensemble Kinetic Energy');
xlabel('Time step');
ylabel('Ensemble kinetic energy');
savefig('EnsembleKETimeSeries.fig')
equlibriumTimeStep = 3000;

% equlibriumTimeStep = ...
%     input('please enter the time step of system equlibrium state:\n');

%% Ensemble Pressure 
% Virial theory : $<\Sigma_{i} f_{i} x_{i}> -2PV = -N_{f} k_{B} T$
ensembleTemperature = 1;

virialFunction = zeros(length(equlibriumTimeStep:M),1);

for i = equlibriumTimeStep:M
   for j = 1:NFX
      virialFunction(i) = FX(i,j)*X(i,j) + FY(i,j)*Y(i,j);
   end
end

ensemblePressure = ...
    (mean(virialFunction) + NFX*problemDim*ensembleTemperature) / problemDim;


%% Particle trace and Ponicare plane
%--------------------------------------------------------------------------------------------------

% Label name cell array
nameChar = cell(1,NVX);
for i  = 1:NVX
      nameChar{i} = ['No.',num2str(i)];
end

%%%% Final let's look how particles' traces like
f = figure('Color',[1 1 1]);
set(f,'Name','Particle traces');
subplot(1,2,1)
plot(X(:,1),Y(:,1),'r','linewidth',2);
hold on;
plot(X(:,2),Y(:,2),'g','linewidth',2);
plot(X(:,3),Y(:,3),'m','linewidth',2);
plot(X(:,4),Y(:,4),'y','linewidth',2);
plot(X(:,5),Y(:,5),'b','linewidth',2);
hold off;
xlabel('Coordinate X');
ylabel('Cooridinate Y');
title('Particles No.1 - No.5 traces');

subplot(1,2,2);
plot(X(:,6),Y(:,6),'r','linewidth',2);
hold on;
plot(X(:,7),Y(:,7),'g','linewidth',2);
plot(X(:,8),Y(:,8),'m','linewidth',2);
plot(X(:,9),Y(:,9),'y','linewidth',2);
plot(X(:,10),Y(:,10),'b','linewidth',2);
hold off;
xlabel('Coordinate X');
ylabel('Cooridinate Y');
title('Particles No.6 - No.10 traces');
savefig('particlesTraces.fig')

%%%% Next we plot the Ponicare plane
% Plot Phase space (Here we decompose the $3\times 4D = 12D$ Phase space into six 2D space)
f = figure('Color',[1 1 1]);
set(f,'Name','Ponicare planes Paticles No.1 - No.5');
% $(\p_{1,x},q_{1,x})$ 
subplot(5,2,1);
plot(X(:,1),VX(:,1),'r')
title('Particle No.1:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{1,y},q_{1,y})$ 
subplot(5,2,2);
plot(Y(:,1),VY(:,1),'r')
title('Particle No.1:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{2,x},q_{2,x})$ 
subplot(5,2,3);
plot(X(:,2),VX(:,2),'b')
title('Particle No.2:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{2,y},q_{2,y})$ 
subplot(5,2,4);
plot(Y(:,2),VY(:,2),'b')
title('Particle No.2:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{3,x},q_{3,x})$ 
subplot(5,2,5);
plot(X(:,3),VX(:,3),'m')
title('Particle No.3:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{3,y},q_{3,y})$ 
subplot(5,2,6);
plot(Y(:,3),VY(:,3),'m')
title('Particle No.3:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{4,x},q_{4,x})$ 
subplot(5,2,7);
plot(X(:,4),VX(:,4),'m')
title('Particle No.4:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{4,y},q_{4,y})$ 
subplot(5,2,8);
plot(Y(:,4),VY(:,4),'m')
title('Particle No.4:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{3,x},q_{3,x})$ 
subplot(5,2,9);
plot(X(:,5),VX(:,5),'m')
title('Particle No.5:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{5,y},q_{5,y})$ 
subplot(5,2,10);
plot(Y(:,5),VY(:,5),'m')
title('Particle No.5:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
savefig('ponicarePlnae1.fig')

f = figure('Color',[1 1 1]);
set(f,'Name','Ponicare planes Paticles No.6 - No.10');
% $(\p_{1,x},q_{1,x})$ 
subplot(5,2,1);
plot(X(:,6),VX(:,6),'r')
title('Particle No.6:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{1,y},q_{1,y})$ 
subplot(5,2,2);
plot(Y(:,6),VY(:,6),'r')
title('Particle No.6:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{2,x},q_{2,x})$ 
subplot(5,2,3);
plot(X(:,7),VX(:,7),'b')
title('Particle No.7:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{2,y},q_{2,y})$ 
subplot(5,2,4);
plot(Y(:,7),VY(:,7),'b')
title('Particle No.7:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{3,x},q_{3,x})$ 
subplot(5,2,5);
plot(X(:,8),VX(:,8),'m')
title('Particle No.8:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{3,y},q_{3,y})$ 
subplot(5,2,6);
plot(Y(:,8),VY(:,8),'m')
title('Particle No.8:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{4,x},q_{4,x})$ 
subplot(5,2,7);
plot(X(:,9),VX(:,9),'m')
title('Particle No.9:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{4,y},q_{4,y})$ 
subplot(5,2,8);
plot(Y(:,9),VY(:,9),'m')
title('Particle No.9:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{3,x},q_{3,x})$ 
subplot(5,2,9);
plot(X(:,10),VX(:,10),'m')
title('Particle No.10:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{5,y},q_{5,y})$ 
subplot(5,2,10);
plot(Y(:,10),VY(:,10),'m')
title('Particle No.10:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
savefig('ponicarePlnae2.fig')


%--------------------------------------------------------------------------------------------------

%% Statistical Analysis and plots

%%%% First let's look at our universe by using boxplot visualization technique
%--------------------------------------------------------------------------------------------------
%%% this part is for velocity
% Velocity X component box plot
f = figure('Color',[1 1 1]);
set(f,'Name','Velocity X components box plot');
boxplot(VX(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('V_{x}');
title('Boxplot of velocity X component box plot');
savefig('boxPlotVx.fig')

% Velocity Y component box plot 
f = figure('Color',[1 1 1]);
set(f,'Name','Velocity Y components box plot');
boxplot(VY(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('V_{y}');
title('Boxplot of velocity Y component');
savefig('boxPlotVy.fig')

%%% this part is for position
% Position X component box plot
f = figure('Color',[1 1 1]);
set(f,'Name','Position X components box plot');
boxplot(X(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('X');
title('Boxplot of position X component');
savefig('boxPlotX.fig')

% Position Y component box plot
f = figure('Color',[1 1 1]);
set(f,'Name','Position Y components box plot');
boxplot(Y(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('Y');
title('Boxplot of position Y component');
savefig('boxPlotY.fig')


%%% this part is for energy
% Kinetic energy
f = figure('Color',[1 1 1]);
set(f,'Name','Box plot of Kinetic energy');
boxplot(Kinetic(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('Kinetic energy');
title('Boxplot of Kinetic energy');
savefig('boxPlotKE.fig')

% Potential energy
f = figure('Color',[1 1 1]);
set(f,'Name','Box plot of Potential energy');
boxplot(Potential(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('Potential energy');
title('Boxplot of Potential energy');
savefig('boxPlotPE.fig')

% Total energy 
f = figure('Color',[1 1 1]);
set(f,'Name','Box plot of Mechainical energy');
boxplot(Energy(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('Mechainical energy');
title('Boxplot of Mechainical energy');
savefig('boxPlotME.fig')

%---------------------------------------------------------------------------------------------------

%%%% Next let's answer questions by using statistics toolbox
%---------------------------------------------------------------------------------------------------
% matrixplot of Velocity of each particle
f = figure('Color',[1 1 1]);
set(f,'Name','Matrixplot of velocity of each particle');
temp = [VX(equlibriumTimeStep:M,1),VY(equlibriumTimeStep:M,1),...
        VX(equlibriumTimeStep:M,2),VY(equlibriumTimeStep:M,2),...
        VX(equlibriumTimeStep:M,3),VY(equlibriumTimeStep:M,3),...
        VX(equlibriumTimeStep:M,4),VY(equlibriumTimeStep:M,4),...
        VX(equlibriumTimeStep:M,5),VY(equlibriumTimeStep:M,5),...
        VX(equlibriumTimeStep:M,6),VY(equlibriumTimeStep:M,6),...
        VX(equlibriumTimeStep:M,7),VY(equlibriumTimeStep:M,7),...
        VX(equlibriumTimeStep:M,8),VY(equlibriumTimeStep:M,8),...
        VX(equlibriumTimeStep:M,9),VY(equlibriumTimeStep:M,9),...
        VX(equlibriumTimeStep:M,10),VY(equlibriumTimeStep:M,10)];    
plotmatrix(temp);
title('Matrixplot of velocity of each particle');
savefig('MatrixplotVelocity.fig');

% matrixplot of potential energy
f = figure('Color',[1 1 1]);
set(f,'Name','Matrixplot of potential energy of each particle');
temp = [Potential(equlibriumTimeStep:M,1),...
        Potential(equlibriumTimeStep:M,2),...
        Potential(equlibriumTimeStep:M,3),...
        Potential(equlibriumTimeStep:M,4),...
        Potential(equlibriumTimeStep:M,5),...
        Potential(equlibriumTimeStep:M,6),...
        Potential(equlibriumTimeStep:M,7),...
        Potential(equlibriumTimeStep:M,8),...
        Potential(equlibriumTimeStep:M,9),...
        Potential(equlibriumTimeStep:M,10)]; 
plotmatrix(temp);
title('Matrixplot of potential energy of each particle');
savefig('MatrixplotPEs.fig');

% boxplot for particle No.1 
f = figure('Color',[1 1 1]);
set(f,'Name','Matrixplot of velocity and position for particle No.1');
temp = [VX(equlibriumTimeStep:M,1),VY(equlibriumTimeStep:M,1),...
         X(equlibriumTimeStep:M,2), Y(equlibriumTimeStep:M,2)];
plotmatrix(temp);
title('Matrixplot of velocity and position for particle No.1');
savefig('MatrixplotVRNo1.fig');

% boxplot for energy
f = figure('Color',[1 1 1]);
set(f,'Name','Matrixplot of ensemble energy');
temp = [Ensemble_Kinetic(equlibriumTimeStep:M),...
        Ensemble_Potential(equlibriumTimeStep:M),...
        Ensemble_Energy(equlibriumTimeStep:M)];
plotmatrix(temp);
title('Matrixplot of ensemble energy');
savefig('MatrixplotEnsembleEnergy.fig');
clear temp

% Maxwell-Boltzmann distrubution verify
% we only research particle 1
temp = VX(equlibriumTimeStep:M,1);
pd = fitdist(temp,'Normal');
x_values = min(temp):(max(temp)-min(temp))/1000:max(temp);
pdfFit = pdf(pd,x_values);
[f_ecdf,x_ecdf] = ecdf(temp);

f = figure('Color',[1 1 1]);
set(f,'Name','Empirical cumulative PDF and PDF fitting');
ecdfhist(f_ecdf,x_ecdf,100);
hold on;
plot(x_values,pdfFit,'r','LineWidth',2);
hold off;
title('Empirical cumulative PDF and PDF fitting');
xlabel('Stochastic variable');
ylabel('Frequence or Probability density');
str1 = ['\mu = ',num2str(pd.mu)];
str2 = ['\sigma = ',num2str(pd.sigma)];
annotation(f,'textbox',...
    [0.174214285714286 0.713028169014085 0.157663443543356 0.106019450033535],...
    'String',[str1,sprintf('\n'),str2],'FitBoxToText','off');
savefig('MaxwellBoltzmannNo1.fig');
clear temp;

% virial function statistics
[f_ecdf,x_ecdf] = ecdf(virialFunction(equlibriumTimeStep:M));
f = figure('Color',[1 1 1]);
set(f,'Name','ECDF density of Virial function');
ecdfhist(f_ecdf,x_ecdf,100);
title('ECDF density of Virial function');
xlabel('Virial function time samples');
ylabel('Frequency');
savefig('virialFunECDFHist.fig');





















%% Program End