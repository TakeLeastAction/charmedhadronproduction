d=xlsread('LHCFlow.xlsx','phi4','L12:N18');  % phi
pt = d(:,1);fp = d(:,2);err=d(:,3);
mt=sqrt(pt.^2+1.02^2);
d=xlsread('LHCFlow.xlsx','Omega4','L12:N19');  % Omega
pt2 = d(:,1);fp2 = d(:,2);err2=d(:,3);
mt2=sqrt(pt2.^2+1.67^2);
figure;subplot(2,1,1);
errorbare('v',mt/2,fp/2,err/2) ;hold on; errorbare('v',mt2/3,fp2/3,err2/3,'r') ;
xlabel('mt/n (GeV/c)');ylabel('v2/n');
axis([0 2 0 0.12]);subplot(2,1,2);
errorbare('v',pt/2,fp/2,err/2) ;hold on; errorbare('v',pt2/3,fp2/3,err2/3,'r') ;
axis([0 2 0 0.12]);xlabel('pt/n (GeV/c)');ylabel('v2/n');
