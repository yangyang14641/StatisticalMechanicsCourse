%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Program name : MDDataAnalysis.m
%%%% Program Author: Yang Yang
%%%% Data : 30/4/2016
%%%% Version: V1.0
%%%% copyright: Author only!
%%%% input file format: result.txt    
%%%% data column ( 
%%%% t, x1, y1, z1, u1, v1, w1, fx1, fy1, fz1, k1, p1, 
%%%%    x2, y2, z2, u2, v2, w2, fx2, fy2, fz2, k2, p2, 
%%%% ......
%%%%    x10, y10, z10, u10, v10, w10, fx10, fy10, fz10, k10, p10)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data input module
filename = 'result.txt';        % data file name goes here
Data = importdata(filename);    % import command
[M,N] = size(Data);             % find dimension of data

problemDim = 3;                 % problem dimension

%% Regroup data
col_num = 11;
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
Y = zeros(M,NY);
for j = 1:NY
    Y(:,j) = Data(:,index + (j-1)*col_num);
end
index = index + 1;

% z-coordinate arrays
NZ = (N-1)/col_num;
Z = zeros(M,NY);
for j = 1:NZ
    Z(:,j) = Data(:,index + (j-1)*col_num);
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

% Velocity z component arrays
NVZ = (N-1)/col_num;
VZ = zeros(M,NVZ);
for j = 1:NVZ
    VZ(:,j) = Data(:,index + (j-1)*col_num);
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

% internal force z component arrays
NFZ = (N-1)/col_num;
FZ = zeros(M,NVZ);
for j = 1:NVZ
    FZ(:,j) = Data(:,index + (j-1)*col_num);
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
plot(Ensemble_Energy);
title('Ensemble total Energy');
xlabel('Time step');
ylabel('Ensemble total energy');
savefig('EnsembleTETimeSeries.fig')
equlibriumTimeStep = 10;

% equlibriumTimeStep = ...
%     input('please enter the time step of system equlibrium state:\n');

%% Ensemble Pressure 
% Virial theory : $<\Sigma_{i} f_{i} x_{i}> -2PV = -N_{f} k_{B} T$

virialFunction = zeros(length(equlibriumTimeStep:M),1);

for i = equlibriumTimeStep:M
   for j = 1:NFX
      virialFunction(i) = FX(i,j)*X(i,j) + FY(i,j)*Y(i,j);
   end
end

ensemblePressure = ...
    (mean(virialFunction) + mean(Ensemble_Kinetic)) / problemDim;

%% Particle intial position
f = figure('Color',[1 1 1]);
set(f,'Name','Particle initial position','Color',[1 1 1]);
[XX,YY,ZZ] = sphere;
surf(XX,YY,ZZ,'FaceAlpha',0.2,'EdgeAlpha',0.1);
hold on;
plot3(X(1,:),Y(1,:),Z(1,:),'r.','MarkerSize',30);
for i = 1:NX
    line([0,X(1,i)],[0,Y(1,i)],[0,Z(1,i)],'linewidth',2);
end
box on;
grid off;
axis equal;
title('Particle initial position');
hold off;

clear XX YY ZZ 


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
plot3(X(:,1),Y(:,1),Z(:,1),'r','linewidth',2);
hold on;
plot3(X(:,2),Y(:,2),Z(:,2),'g','linewidth',2);
plot3(X(:,3),Y(:,3),Z(:,3),'m','linewidth',2);
plot3(X(:,4),Y(:,4),Z(:,4),'y','linewidth',2);
plot3(X(:,5),Y(:,5),Z(:,5),'b','linewidth',2);
hold off;
box on;
xlabel('Coordinate X');
ylabel('Cooridinate Y');
zlabel('Cooridinate Z');
title('Particles No.1 - No.5 traces');

subplot(1,2,2);
plot3(X(:,6),Y(:,6),Z(:,6),'r','linewidth',2);
hold on;
plot3(X(:,7),Y(:,7),Z(:,7),'g','linewidth',2);
plot3(X(:,8),Y(:,8),Z(:,8),'m','linewidth',2);
plot3(X(:,9),Y(:,9),Z(:,9),'y','linewidth',2);
plot3(X(:,10),Y(:,10),Z(:,10),'b','linewidth',2);
hold off;
box on;
xlabel('Coordinate X');
ylabel('Cooridinate Y');
zlabel('Cooridinate Z');
title('Particles No.6 - No.10 traces');
savefig('particlesTraces.fig')

%%%% Next we plot the Ponicare plane
% Plot Phase space (Here we decompose the $3\times 4D = 12D$ Phase space into six 2D space)
f = figure('Color',[1 1 1]);
set(f,'Name','Ponicare planes Paticles No.1 - No.5');
% $(\p_{1,x},q_{1,x})$ 
subplot(5,3,1);
plot(X(:,1),VX(:,1),'r')
title('Particle No.1:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{1,y},q_{1,y})$ 
subplot(5,3,2);
plot(Y(:,1),VY(:,1),'r')
title('Particle No.1:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{1,z},q_{1,z})$ 
subplot(5,3,3);
plot(Z(:,1),VZ(:,1),'r')
title('Particle No.1:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');

% $(\p_{2,x},q_{2,x})$ 
subplot(5,3,4);
plot(X(:,2),VX(:,2),'r')
title('Particle No.2:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{2,y},q_{2,y})$ 
subplot(5,3,5);
plot(Y(:,2),VY(:,2),'r')
title('Particle No.2:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{2,z},q_{2,z})$ 
subplot(5,3,6);
plot(Z(:,2),VZ(:,2),'r')
title('Particle No.2:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');

% $(\p_{3,x},q_{3,x})$ 
subplot(5,3,7);
plot(X(:,3),VX(:,3),'r')
title('Particle No.3:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{3,y},q_{3,y})$ 
subplot(5,3,8);
plot(Y(:,3),VY(:,3),'r')
title('Particle No.3:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{3,z},q_{3,z})$ 
subplot(5,3,9);
plot(Z(:,3),VZ(:,3),'r')
title('Particle No.3:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');

% $(\p_{4,x},q_{4,x})$ 
subplot(5,3,10);
plot(X(:,4),VX(:,4),'r')
title('Particle No.4:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{4,y},q_{4,y})$ 
subplot(5,3,11);
plot(Y(:,4),VY(:,4),'r')
title('Particle No.4:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{4,z},q_{4,z})$ 
subplot(5,3,12);
plot(Z(:,4),VZ(:,4),'r')
title('Particle No.4:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');

% $(\p_{5,x},q_{5,x})$ 
subplot(5,3,13);
plot(X(:,5),VX(:,5),'r')
title('Particle No.5:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{5,y},q_{5,y})$ 
subplot(5,3,14);
plot(Y(:,5),VY(:,5),'r')
title('Particle No.5:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{5,z},q_{5,z})$ 
subplot(5,3,15);
plot(Z(:,5),VZ(:,5),'r')
title('Particle No.5:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');
savefig('ponicarePlnae1.fig')


f = figure('Color',[1 1 1]);
set(f,'Name','Ponicare planes Paticles No.6 - No.10');
% $(\p_{6,x},q_{6,x})$ 
subplot(5,3,1);
plot(X(:,6),VX(:,6),'r')
title('Particle No.6:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{6,y},q_{6,y})$ 
subplot(5,3,2);
plot(Y(:,6),VY(:,6),'r')
title('Particle No.6:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{6,z},q_{6,z})$ 
subplot(5,3,3);
plot(Z(:,6),VZ(:,6),'r')
title('Particle No.6:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');

% $(\p_{7,x},q_{7,x})$ 
subplot(5,3,4);
plot(X(:,7),VX(:,7),'r')
title('Particle No.7:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{7,y},q_{7,y})$ 
subplot(5,3,5);
plot(Y(:,7),VY(:,7),'r')
title('Particle No.7:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{7,z},q_{7,z})$ 
subplot(5,3,6);
plot(Z(:,7),VZ(:,7),'r')
title('Particle No.7:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');

% $(\p_{8,x},q_{8,x})$ 
subplot(5,3,7);
plot(X(:,8),VX(:,8),'r')
title('Particle No.8:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{8,y},q_{8,y})$ 
subplot(5,3,8);
plot(Y(:,8),VY(:,8),'r')
title('Particle No.8:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{8,z},q_{8,z})$ 
subplot(5,3,9);
plot(Z(:,8),VZ(:,8),'r')
title('Particle No.8:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');

% $(\p_{9,x},q_{9,x})$ 
subplot(5,3,10);
plot(X(:,9),VX(:,9),'r')
title('Particle No.9:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{9,y},q_{9,y})$ 
subplot(5,3,11);
plot(Y(:,9),VY(:,9),'r')
title('Particle No.9:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{9,z},q_{9,z})$ 
subplot(5,3,12);
plot(Z(:,9),VZ(:,9),'r')
title('Particle No.9:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');

% $(\p_{10,x},q_{10,x})$ 
subplot(5,3,13);
plot(X(:,10),VX(:,10),'r')
title('Particle No.10:  q_x V.s p_x');
xlabel('q_x');
ylabel('p_x');
% $(\p_{10,y},q_{10,y})$ 
subplot(5,3,14);
plot(Y(:,10),VY(:,10),'r')
title('Particle No.10:  q_y V.s p_y');
xlabel('q_y');
ylabel('p_y');
% $(\p_{10,z},q_{10,z})$ 
subplot(5,3,15);
plot(Z(:,10),VZ(:,10),'r')
title('Particle No.10:  q_z V.s p_z');
xlabel('q_z');
ylabel('p_z');
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

% Velocity Z component box plot
f = figure('Color',[1 1 1]);
set(f,'Name','Velocity Z components box plot');
boxplot(VZ(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('V_{z}');
title('Boxplot of velocity Z component box plot');
savefig('boxPlotVz.fig')

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

% Position Z component box plot
f = figure('Color',[1 1 1]);
set(f,'Name','Position Z components box plot');
boxplot(Z(equlibriumTimeStep:M,:),'labels',nameChar);
ylabel('Z');
title('Boxplot of position Z component');
savefig('boxPlotZ.fig')


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
% f = figure('Color',[1 1 1]);
% set(f,'Name','Matrixplot of velocity of each particle');
% temp = [VX(equlibriumTimeStep:M,1),VY(equlibriumTimeStep:M,1),VY(equlibriumTimeStep:M,1),...
%         VX(equlibriumTimeStep:M,2),VY(equlibriumTimeStep:M,2),VY(equlibriumTimeStep:M,2),...
%         VX(equlibriumTimeStep:M,3),VY(equlibriumTimeStep:M,3),VY(equlibriumTimeStep:M,3),...
%         VX(equlibriumTimeStep:M,4),VY(equlibriumTimeStep:M,4),VY(equlibriumTimeStep:M,4),...
%         VX(equlibriumTimeStep:M,5),VY(equlibriumTimeStep:M,5),VY(equlibriumTimeStep:M,5),...
%         VX(equlibriumTimeStep:M,6),VY(equlibriumTimeStep:M,6),VY(equlibriumTimeStep:M,6),...
%         VX(equlibriumTimeStep:M,7),VY(equlibriumTimeStep:M,7),VY(equlibriumTimeStep:M,7),...
%         VX(equlibriumTimeStep:M,8),VY(equlibriumTimeStep:M,8),VY(equlibriumTimeStep:M,8),...
%         VX(equlibriumTimeStep:M,9),VY(equlibriumTimeStep:M,9),VY(equlibriumTimeStep:M,9),...
%         VX(equlibriumTimeStep:M,10),VY(equlibriumTimeStep:M,10),VY(equlibriumTimeStep:M,10)];    
% plotmatrix(temp);
% title('Matrixplot of velocity of each particle');
% savefig('MatrixplotVelocity.fig');

% matrixplot of Kinetic energy
f = figure('Color',[1 1 1]);
set(f,'Name','Matrixplot of potential energy of each particle');
temp = [Kinetic(equlibriumTimeStep:M,1),...
        Kinetic(equlibriumTimeStep:M,2),...
        Kinetic(equlibriumTimeStep:M,3),...
        Kinetic(equlibriumTimeStep:M,4),...
        Kinetic(equlibriumTimeStep:M,5),...
        Kinetic(equlibriumTimeStep:M,6),...
        Kinetic(equlibriumTimeStep:M,7),...
        Kinetic(equlibriumTimeStep:M,8),...
        Kinetic(equlibriumTimeStep:M,9),...
        Kinetic(equlibriumTimeStep:M,10)]; 
plotmatrix(temp);
title('Matrixplot of kinetic energy of each particle');
savefig('MatrixplotKEs.fig');

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
temp = [VX(equlibriumTimeStep:M,1),VY(equlibriumTimeStep:M,1),VZ(equlibriumTimeStep:M,1),...
         X(equlibriumTimeStep:M,2), Y(equlibriumTimeStep:M,2),Z(equlibriumTimeStep:M,1)];
plotmatrix(temp);
title('Matrixplot of velocity and position for particle No.1');
savefig('MatrixplotVRNo1.fig');

% matrixplot for energy
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
% temp = VX(equlibriumTimeStep:M,1);
% pd = fitdist(temp,'Normal');
% x_values = min(temp):(max(temp)-min(temp))/1000:max(temp);
% pdfFit = pdf(pd,x_values);
% [f_ecdf,x_ecdf] = ecdf(temp);
% 
% f = figure('Color',[1 1 1]);
% set(f,'Name','Empirical cumulative PDF and PDF fitting');
% ecdfhist(f_ecdf,x_ecdf,100);
% hold on;
% plot(x_values,pdfFit,'r','LineWidth',2);
% hold off;
% title('Empirical cumulative PDF and PDF fitting');
% xlabel('Stochastic variable');
% ylabel('Frequence or Probability density');
% str1 = ['\mu = ',num2str(pd.mu)];
% str2 = ['\sigma = ',num2str(pd.sigma)];
% annotation(f,'textbox',...
%     [0.174214285714286 0.713028169014085 0.157663443543356 0.106019450033535],...
%     'String',[str1,sprintf('\n'),str2],'FitBoxToText','off');
% savefig('MaxwellBoltzmannNo1.fig');
% clear temp;
% 

% virial function statistics
[f_ecdf,x_ecdf] = ecdf(virialFunction(equlibriumTimeStep:M));
f = figure('Color',[1 1 1]);
set(f,'Name','ECDF density of Virial function');
ecdfhist(f_ecdf,x_ecdf,100);
xlabel('Virial function time samples');
ylabel('Frequency');
title('ECDF density of Virial function');
savefig('virialFunECDFHist.fig');

















%% Program End
