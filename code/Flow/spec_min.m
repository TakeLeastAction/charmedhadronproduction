function [x,fval,exitflag,output,lambda,grad,hessian] = spec_min(x0,aineq,bineq,fGoal)
% This is an auto generated M-file from Optimization Tool.

% Start with the default options
options = optimset;
% Modify options setting
options = optimset(options,'Display', 'off');
options = optimset(options,'Algorithm', 'active-set');

%options = optimset(options,'TolFun', 1E-3);

%options = optimset(options,'TolX', 1E-4);

[x,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon(@(x)fGoal(x),x0,aineq,bineq,[],[],[],[],[],options);
