function  [ best_T,yield, av_PT, vel]=scan_noas_opt(Data1,m0)
global   Data  m 
Data=Data1;
m=m0;
pt=Data(1,:);fp=Data(2,:);err=Data(3,:);

x0=[0.08    10.5    0.8].'; % T yidu  beta0
%x0=[1.178500e-01 3.992912e+00 9.591875e-01];


[x,fval,exitflag,output,lambda,grad,hessian] = spec_min(x0,Aineq,Bineq,@Goalf);
[x1,fval1,exitflag,output]=patternsearch(@Goalf,x0,Aineq,Bineq);
%%
best_T=x(1);
best_xi=x(2);
best_beta=x(3);
%% best fit
SPT=Func1(m0,pt,best_T,best_beta,best_xi);
SPT=SPT.';
plot(pt,SPT,'-r.')
yield=trapz(pt,pt.*SPT)*2*pi;
av_PT=trapz(pt,pt.*SPT.*pt)*2*pi/yield;
vel=av_PT/sqrt(m^2+av_PT^2);


