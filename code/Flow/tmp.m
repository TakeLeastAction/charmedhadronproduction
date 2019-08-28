
pt=linspace(0,6,20);

 
xi_s=0.785;
xi_anti_s=xi_s;      %/exp(2*10/175);
xi_c=1.0;
xi_anti_c=1.0;

HOSTNAME='bl-1-1.physics.sjtu.edu.cn';
USERNAME='sunkj';
PASSWORD='N123456';
ssh2_conn = ssh2_config(HOSTNAME,USERNAME,PASSWORD);
%input=[0.175,beta,R0,tau0, dtau, m_u, m_s, m_c, 1.0, 1.0, 1.0, 1.0,1.0, 1.0, 0.0 1.0 0.0 0.0 0.0 0.0 0.0].';
ell_1=0.0;ell_2=0.0; ell_3=0.0;s4=0.0;
%input=[0.14,beta,R0,tau0, dtau, m_u, m_s, m_c, 1.0, 1.0, 1.0, 1.0,1.0, 1.0, c1, c2, s2, ell_1, ell_2, ell_3, s4].';
%save input.txt -ascii input

%%
file1='result_phi.dat';
file2='result_Omega.dat';
%file3='result_p.txt';
%file4='result_lambda.txt';
%file5='result_xi.txt';
%file6='result_Occc.txt';
local_path1=pwd;
remote_path1='/home/sunkj/Fitting/tmpLHCcharmflow';
ssh4 = scp_get(ssh2_conn,{file1},local_path1,remote_path1);
SPT_phi=importdata('result_phi.dat','',0);
SPT_Omega=importdata('result_Omega.dat','',0);



Nf=25;
N= length(pt);
faip=linspace(0,pi*2,25).';  % need to be the same  in the cluster coalsescence code.
%mt=sqrt(pt.^2+4*m0^2);
for k=1:N
    fp_phi(k)=trapz(faip,SPT_phi((k-1)*Nf+1:k*Nf));
    v2_phi(k)=trapz(faip,SPT_phi((k-1)*Nf+1:k*Nf).*cos(2*faip))./trapz(faip,SPT_phi((k-1)*Nf+1:k*Nf));
    
    %fp_Omega(k)=trapz(faip,SPT_Omega((k-1)*Nf+1:k*Nf));
    %v2_Omega(k)=trapz(faip,SPT_Omega((k-1)*Nf+1:k*Nf).*cos(2*faip))./trapz(faip,SPT_Omega((k-1)*Nf+1:k*Nf));
    
end
SPT_phi=pt.'.*fp_phi*3/4*xi_s*xi_anti_s;
%SPT_Omega=pt.'.*fp_Omega*4/8*xi_s^3;

figure;plot(pt,v2_phi,'-r.')