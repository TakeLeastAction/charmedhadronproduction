function [SPT_phi,v2_phi,SPT_Omega,v2_Omega]=Func3(pt,param)
%function [SPT_phi,SPT_Omega,SPT_p,SPT_lambda,SPT_xi]=Func3(pt,param)
% dN/2*pi*ptdptd(yita)
%cmd='cd /home/sunkj/compare/matlab_serv;ls;rm input.txt;ls ';
%ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD);
%ssh2_conn = ssh2_command(ssh2_conn, cmd);
% or to suppress the output from being displayed I can run
% Reusing Advanced Connection:
%ssh2_conn = ssh2_command(ssh2_conn, 'ls -la *compare*',1);
%command_response = ssh2_command_response(ssh2_conn);
%global SPT_Occc

save PT.txt -ascii pt;
fid = fopen('jobNum.txt','wt');
fprintf(fid,'%g\n',length(pt));
fclose(fid);


beta=param(1);
tau0=param(2);
dtau = param(3);
R0 = 7+0.6*tau0;
m_u=0.3;
m_s=0.5;
m_c=1.5;
xi_u=param(5);
xi_anti_u=xi_u ;   %/exp(2*20/3/175);
xi_s=param(4);
xi_anti_s=xi_s;      %/exp(2*10/175);
xi_c=1.0;
xi_anti_c=1.0;
c1 = param(6);
c2 = param(7);
s2 = param(8);
HOSTNAME='202.120.47.121';
USERNAME='pujie';
PASSWORD='N123456';
ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD);
%input=[0.175,beta,R0,tau0, dtau, m_u, m_s, m_c, 1.0, 1.0, 1.0, 1.0,1.0, 1.0, 0.0 1.0 0.0 0.0 0.0 0.0 0.0].';
ell_1=0.0;ell_2=0.0; ell_3=0.0;s4=0.0;
input=[0.14,beta,R0,tau0, dtau, m_u, m_s, m_c, 1.0, 1.0, 1.0, 1.0,1.0, 1.0, c1, c2, s2, ell_1, ell_2, ell_3, s4].';
save input.txt -ascii input

file_upload='input.txt';
local_path0=pwd;
remote_path0='/home/pujie/Fitting/LHCFlow';
cmd='cd /home/pujie/Fitting/LHCFlow;rm input.txt result*.txt';
% remove input.txt on server
ssh1 = ssh2_command(ssh2_conn,cmd);
% upload input.txt on server
ssh2 = scp_put(ssh1,file_upload,remote_path0,local_path0);
remote_path2='/home/pujie/Fitting/LHCFlow/1phi';
remote_path3='/home/pujie/Fitting/LHCFlow/2Omega';
%remote_path4='/home/sunkj/Fitting/LHC/3proton';
%remote_path5='/home/sunkj/Fitting/LHC/4Lambda';
%remote_path6='/home/sunkj/Fitting/LHC/5Xi';
%remote_path7='/home/sunkj/Fitting/LHC/6Occc';

scp_put(ssh1,{'PT.txt','jobNum.txt'},remote_path2,local_path0)
scp_put(ssh1,{'PT.txt','jobNum.txt'},remote_path3,local_path0)
%scp_put(ssh1,{'PT.txt','jobNum.txt'},remote_path4,local_path0)
%scp_put(ssh1,{'PT.txt','jobNum.txt'},remote_path5,local_path0)
%scp_put(ssh1,{'PT.txt','jobNum.txt'},remote_path6,local_path0)
%scp_put(ssh1,{'PT.txt','jobNum.txt'},remote_path7,local_path0)

%%
% run python_serv.py
%cmd='cd /home/sunkj/Fitting/LHCFlow;python Coal_submit2.py';
%ssh3 = ssh2_command(ssh2_conn,cmd);

%%
file1='result_phi.dat';
file2='result_Omega.dat';
%file3='result_p.txt';
%file4='result_lambda.txt';
%file5='result_xi.txt';
%file6='result_Occc.txt';
local_path1=pwd;
remote_path1='/home/pujie/Fitting/LHCFlow';
ssh4 = scp_get(ssh1,{file1,file2},local_path1,remote_path1);
SPT_phi=importdata('result_phi.dat','',0);
SPT_Omega=importdata('result_Omega.dat','',0);

Nf=25;
N= length(pt);
faip=linspace(0,pi*2,25).';  % need to be the same  in the cluster coalsescence code.
%mt=sqrt(pt.^2+4*m0^2);
for k=1:N
    fp_phi(k)=trapz(faip,SPT_phi((k-1)*Nf+1:k*Nf));
    v2_phi(k)=trapz(faip,SPT_phi((k-1)*Nf+1:k*Nf).*cos(2*faip))./trapz(faip,SPT_phi((k-1)*Nf+1:k*Nf));
    
    fp_Omega(k)=trapz(faip,SPT_Omega((k-1)*Nf+1:k*Nf));
    v2_Omega(k)=trapz(faip,SPT_Omega((k-1)*Nf+1:k*Nf).*cos(2*faip))./trapz(faip,SPT_Omega((k-1)*Nf+1:k*Nf));
    
end
SPT_phi=pt.'.*fp_phi*3/4*xi_s*xi_anti_s;
SPT_Omega=pt.'.*fp_Omega*4/8*xi_s^3;
%SPT_p=importdata('result_p.txt','',0);
%SPT_p=2*pi*pt.*SPT_p*2/8*xi_u^3;
%SPT_lambda=importdata('result_lambda.txt','',0);
%SPT_lambda=2*pi*pt.*SPT_lambda*2/8*xi_u^2*xi_s;
%SPT_xi=importdata('result_xi.txt','',0);
%SPT_xi=2*pi*pt.*SPT_xi*2/8*xi_u*xi_s^2;
%SPT_Occc=importdata('result_Occc.txt','',0);
%SPT_Occc=2*pi*pt.*SPT_Occc*4/8*xi_c^3;

ssh2_close(ssh2_conn);
ssh2_close(ssh1);
ssh2_close(ssh2);
%ssh2_close(ssh3);
ssh2_close(ssh4);
