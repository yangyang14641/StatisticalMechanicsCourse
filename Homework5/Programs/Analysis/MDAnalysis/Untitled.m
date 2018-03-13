n = 1000;
rho = -.8;
nu = 5;
T = mvtrnd([1 rho; rho 1], nu, n);
U = tcdf(T,nu);
X = [invCDF1(ceil(n1*U(:,1))) invCDF2(ceil(n2*U(:,2)))];

[n1,ctr1] = hist(X(:,1),10);
[n2,ctr2] = hist(X(:,2),10);
subplot(2,2,2);
plot(X(:,1),X(:,2),'.');
axis([-3.5 3.5 -3.5 3.5]);
h1 = gca;
title('1000 Simulated Dependent Values');
xlabel('X1');
ylabel('X2');
subplot(2,2,4);
bar(ctr1,-n1,1);
axis([-3.5 3.5 -max(n1)*1.1 0]);
axis('off');
h2 = gca;
subplot(2,2,1);
barh(ctr2,-n2,1);
axis([-max(n2)*1.1 0 -3.5 3.5]);
axis('off');
h3 = gca;
h1.Position = [0.35 0.35 0.55 0.55];
h2.Position = [.35 .1 .55 .15];
h3.Position = [.1 .35 .15 .55];
colormap([.8 .8 1]);


%% distribution fit function
load hospital
x = hospital.Weight;
pd = fitdist(x,'Normal');

x_values = 50:1:250;
y = pdf(pd,x_values);
plot(x_values,y,'LineWidth',2)
cdfhist()





