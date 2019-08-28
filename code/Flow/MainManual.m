clc;clear;
global   Ept_phi Ev2_phi Eerr_phi  Ept_Omega  Ev2_Omega  Eerr_Omega  Ept_p  Espt_p  Eerr_p  Ept_Lambda  Espt_Lambda Eerr_Lambda  Ept_Xi  Espt_Xi Eerr_Xi
global  tau_max beta_max xi_s_max   tau_min  beta_min xi_s_min  xi_u_max xi_u_min  IT
global SPT_Occc dtau_min dtau_max
tau_max=15.0;tau_min=6.0;beta_max=1.1;beta_min=0.5;xi_s_max=3;xi_s_min=0.5;
xi_u_max=1.0;xi_u_min=1.0; dtau_min =0.5;dtau_max = 1.5;
[Ept_phi, Ev2_phi,Eerr_phi, Ept_Omega, Ev2_Omega,Eerr_Omega]=Initial();  % flow 
save ExpData
IT=1;
tic
%%
c1 = 0.46;%2.2 2.50.38 0.38 0.3 0.2 0.15 0.06 0.07 0.08
c2 = 9.2;%4,0.4 0.7 4 4 9 8 10 11.5 11 10.5 12 10
s2 = -0.025;  %-0.1 -0.15 -0.02 %-0.1%-0.085;
xi_s=1.05;
x0=[1.09,15,1.3,xi_s,1, c1, c2, s2].';%   beta0 tau dtau xi_s xi_u

param=x0;%[1.0,5.0  ,1.3,2.0*1.7].';%[0.795,7.3,1.0,3.0].';
pt=linspace(0,6,20);
pt=pt.';
[SPT_phi,v2_phi,SPT_Omega,v2_Omega]=Func3(pt,param);

% SPT_xi = SPT_xi*3;%% Important!
% SPT_p = SPT_p*5;
% SPT_lambda = SPT_lambda*7.44;

load ExpData
figure;
subplot(2,1,1)
errorbare('v',Ept_phi,Ev2_phi,Eerr_phi,'-c') ;hold on;
hold on;plot(pt,v2_phi,'-r',...
    'LineWidth',3,...
    'Color',[0 0.447058826684952 0.74117648601532]);legend('\phi')

subplot(2,1,2)
errorbare('v',Ept_Omega,Ev2_Omega,Eerr_Omega,'-b') ;hold on;
plot(pt,v2_Omega,'-g',...
    'LineWidth',3,...
    'Color',[0 0.447058826684952 0.74117648601532]);legend('\Omega^-');

figure; errorbare('v',Ept_phi/2,Ev2_phi/2,Eerr_phi/2,'-c') ;hold on;errorbare('v',Ept_Omega/3,Ev2_Omega/3,Eerr_Omega/3,'-b') ;
plot(pt/2,v2_phi/2,'-c');plot(pt/3,v2_Omega/3,'-b');

R0=7+0.6*param(2);
SPT=xi_s*proton_noas(0.5,pt,0.154,param(1),param(2),R0,0.0,c1,c2,s2,0.0,0.0,0.0,0.0,0);
SPT2=xi_s*proton_noas(0.5,pt,0.154,param(1),param(2),R0,0.0,c1,c2,s2,0.0,0.0,0.0,0.0,2);
v2 = SPT2./SPT;
plot(pt,v2,'-r.');

figure; errorbare('v',sqrt(Ept_phi.^2+1.0^2)/2,Ev2_phi/2,Eerr_phi/2,'-c') ;
hold on;errorbare('v',sqrt(Ept_Omega.^2+1.67^2)/3,Ev2_Omega/3,Eerr_Omega/3,'-b') ;
plot(sqrt(pt.^2+1.0^2)/2,v2_phi/2,'-c');plot(sqrt(pt.^2+1.5^2)/3,v2_Omega/3,'-b');plot(sqrt(pt.^2+0.5^2),v2,'-r.');

figure;
subplot(2,1,1)
semilogy(pt,SPT_phi,'-c');legend('\phi');
subplot(2,1,2)
semilogy(pt,SPT_Omega,'-g');legend('\Omega^-');

Nphi = trapz(pt,SPT_phi)
No = trapz(pt,SPT_Omega)