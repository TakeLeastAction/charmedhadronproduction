function  Main
clc;clear;
global   Ept_phi Espt_phi Eerr_phi  Ept_Omega  Espt_Omega  Eerr_Omega  Ept_p  Espt_p  Eerr_p  Ept_Lambda  Espt_Lambda Eerr_Lambda  Ept_Xi  Espt_Xi Eerr_Xi
global  tau_max beta_max xi_s_max   tau_min  beta_min xi_s_min  xi_u_max xi_u_min  IT
global SPT_Occc dtau_min dtau_max
tau_max=10.0;tau_min=8.0;beta_max=1.1;beta_min=0.9;xi_s_max=0.8;xi_s_min=0.5;
xi_u_max=1.0;xi_u_min=0.5; dtau_min =0.5;dtau_max = 1.5;
[Ept_phi, Espt_phi,Eerr_phi, Ept_Omega, Espt_Omega,Eerr_Omega,Ept_p, Espt_p,Eerr_p, Ept_Lambda, Espt_Lambda,Eerr_Lambda, Ept_Xi, Espt_Xi,Eerr_Xi]=Initial();
% figure;errorbare('vlogy',Ept_p,Espt_p,Eerr_p) ;hold on;
% errorbare('vlogy',Ept_Lambda,Espt_Lambda,Eerr_Lambda,'-r') ;errorbare('vlogy',Ept_phi,Espt_phi,Eerr_phi,'-c') ;errorbare('vlogy',Ept_Xi,Espt_Xi,Eerr_Xi,'-k') ;errorbare('vlogy',Ept_Omega,Espt_Omega,Eerr_Omega,'-b') ;
% legend('p','','\Lambda','','\phi','','\Xi','','\Omega^-')

save ExpData

IT=1;
tic
%%
%x0=[0.98,10,1.0,0.6,0.7].';%   beta0 tau dtau xi_s xi_u
x0=[1.007,9,1.055,0.58,0.7469].'; %%   beta0 tau dtau xi_s xi_u
[x,fval,exitflag,output]=patternsearch(@Goalf,x0,[],[],[],[],[beta_min tau_min dtau_min xi_s_min xi_u_min],[beta_max tau_max dtau_max xi_s_max xi_u_max]);
x
fval
% subject=strcat('Hi master: This is the fitting of LHC, the final result  is: ',' objF = ',num2str(fval),', param=',...
%     num2str(x.'));
% content='   Kiss while your lips are red!';
% text=strcat(subject,content);
% report2me(text);


%%
%open Expdata_LHC.fig
param=x;%[1.0,5.0  ,1.3,2.0*1.7].';%[0.795,7.3,1.0,3.0].';
pt=linspace(0,6,30);
pt=pt.';
[SPT_phi,SPT_Omega,SPT_p,SPT_lambda,SPT_xi]=Func3(pt,param);

SPT_xi = SPT_xi*3;%% Important!
SPT_p = SPT_p*5;
SPT_lambda = SPT_lambda*7.44;
load ExpData
figure;
subplot(2,3,1)
errorbare('vlogy',Ept_p,Espt_p,Eerr_p) ;hold on;
plot(pt,SPT_p,'-b');legend('p')
subplot(2,3,2)
errorbare('vlogy',Ept_Lambda,Espt_Lambda,Eerr_Lambda,'-r') ;hold on;
plot(pt,SPT_lambda,'-r');legend('\Lambda')
subplot(2,3,3)
errorbare('vlogy',Ept_phi,Espt_phi,Eerr_phi,'-c') ;hold on;
hold on;plot(pt,SPT_phi,'-c');legend('\phi')
subplot(2,3,4)
errorbare('vlogy',Ept_Xi,Espt_Xi,Eerr_Xi,'-k') ;hold on;
plot(pt,SPT_xi,'-k');legend('\Xi')
subplot(2,3,5)
errorbare('vlogy',Ept_Omega,Espt_Omega,Eerr_Omega,'-b') ;hold on;
plot(pt,SPT_Omega,'-g');legend('\Omega^-');








