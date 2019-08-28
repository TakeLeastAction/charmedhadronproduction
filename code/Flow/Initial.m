function [Ept_phi, Espt_phi,Eerr_phi, Ept_Omega, Espt_Omega,Eerr_Omega]=Initial()

k=1;
while k==1
    try
        %Dp=xlsread('LHCData.xlsx','LHCproton','S8:U42');  %proton  0-5%
        %DL=xlsread('LHCData.xlsx','LHCLambda','K8:M31');  %Lambda 0-5%
        Dphi=xlsread('LHCdatafit.xlsx','1','A4:C10');  %phi 10-20%
%         Dphi=xlsread('LHCdatafit.xlsx','1','A4:C9');  %phi 50-60%
        %DXi=xlsread('LHCData.xlsx','LHCXi','O6:Q23');  %Xi  0-10%
        DO=xlsread('LHCdatafit.xlsx','1','F4:H11');  %Omega^- 10-20%
%         DO=xlsread('LHCdatafit.xlsx','1','F4:H10');  %Omega^- 50-60%
        k=2;
    catch
        disp('reading fails')
        pause(10);
        k=1;
    end
end
%Ept_p = Dp(:,1);Espt_p=Dp(:,2);Eerr_p = Dp(:,3);
Ept_phi = Dphi(:,1);Espt_phi=Dphi(:,2);Eerr_phi = Dphi(:,3);
Ept_Omega = DO(:,1);Espt_Omega=DO(:,2);Eerr_Omega = DO(:,3);  % be careful
%Ept_Lambda = DL(:,1);Espt_Lambda=DL(:,2);Eerr_Lambda = DL(:,3);
%Ept_Xi = DXi(:,1);Espt_Xi=DXi(:,2);Eerr_Xi = DXi(:,3);
% figure;errorbare('vlogy',Ept_p,Espt_p,Eerr_p) ;hold on;
% errorbare('vlogy',Ept_Lambda,Espt_Lambda,Eerr_Lambda,'-r') ;errorbare('vlogy',Ept_phi,Espt_phi,Eerr_phi,'-c') ;errorbare('vlogy',Ept_Xi,Espt_Xi,Eerr_Xi,'-k') ;errorbare('vlogy',Ept_Omega,Espt_Omega,Eerr_Omega,'-b') ;
% legend('p','','\Lambda','','\phi','','\Xi','','(\Omega^-+\Omega^+)/2')
% Ept_p=[0.3240    0.3740    0.4250    0.4750    0.5250    0.5750    0.6260  ...
%     0.6760    0.7260    0.7760    0.8260    0.8770    0.9270    0.9770    ...
%     1.0500    1.1500    1.2500    1.3500    1.4500 1.5500    1.6500    1.7500    ...
%     1.8500    1.9500    2.0500    2.1500    2.2500    2.3500    2.4500    2.5500    2.6500    2.7500].';
% Espt_p=[1.8400    2.0300    2.3100    2.3900    2.7100    2.9900    3.0900 ...
%     3.2900    3.4000    3.6300    3.6300    3.7400    3.7400    3.8700    ...
%     3.8700    3.8700    3.8700    3.6300    3.4000 3.1900    2.8600    2.5400 ...
%     2.2600    1.9900    1.7500    1.5800    1.3300    1.1100    0.9340    0.8100    0.6770    0.5660].';
% 
% Eerr_p=Espt_p/10;
% 
% Ept_phi=[0.6550    0.9070    1.2400    1.7500    2.2600    2.7500    3.5100    4.4500].';
% Espt_phi=[ 6.0500    7.3900    6.6900    4.5800    2.4400    1.0600    0.3100    0.0494].';
% Eerr_phi=Espt_phi/10;
% 
% 
% Ept_Omega=[ 1.3100    1.5000    1.6900    1.9000    2.2000    2.6000    3.0000    3.4000    3.8000    4.2500 ].';
% Espt_Omega=[0.2610    0.2610    0.2490    0.2230    0.1790    0.1170    0.0717    0.0470    0.0240    0.0107].';
% Eerr_Omega=Espt_Omega/10;
% 
% Ept_Lambda=[ 0.6027    0.6964    0.7902    0.8839    0.9911    1.1049    ...
%     1.2188    1.3326    1.4397    1.5402    1.6406    1.7411    1.8616    ...
%     1.9821    2.1295    2.3036    2.4844    2.6987    2.8929 3.0938    3.3013    3.4955    3.6964    3.8973    4.2455].';
% Espt_Lambda=[8.7152    9.6359   10.8265   11.5884   11.9857   12.6116  ...
%     12.8407   12.6116   12.2316   11.6157   11.0684   10.0165    8.8959  ...
%     7.7816    6.5738    5.0422    3.8097    2.5360    1.8496 1.2590    0.8584...
%     0.5944    0.3921    0.2807    0.1539].';
% Eerr_Lambda=Espt_Lambda/10.0;
% 
% Ept_Xi=[0.6820    0.8790    1.0700    1.2600    1.4600    1.6500    1.8400 ...
%     2.0300    2.2300    2.4200    2.6200    2.8100    3.0100    3.2000    3.4000    3.5800    3.7800    3.9700    4.1600].';
% Espt_Xi=[ 3.2500    5.0000    5.2000    5.2600    5.2100    4.4700    ...
%     3.7800    3.3100    2.6100    2.0700    1.5100    1.1300    0.7940    0.5410    0.3720    0.2550    0.1790    0.1330    0.0843].';
% Eerr_Xi=Espt_Xi/10;






cd 1phi
save PT.txt -ascii Ept_phi;
fid = fopen('jobNum.txt','wt'); 
fprintf(fid,'%g\n',length(Ept_phi));       
fclose(fid);
cd ..
cd 2Omega
save PT.txt -ascii Ept_Omega;
fid = fopen('jobNum.txt','wt'); 
fprintf(fid,'%g\n',length(Ept_Omega));       
fclose(fid);
cd ..
% cd 3proton
% save PT.txt -ascii Ept_p;
% fid = fopen('jobNum.txt','wt'); 
% fprintf(fid,'%g\n',length(Ept_p));       
% fclose(fid);
% cd ..
% cd 4Lambda
% save PT.txt -ascii Ept_Lambda;
% fid = fopen('jobNum.txt','wt'); 
% fprintf(fid,'%g\n',length(Ept_Lambda));       
% fclose(fid);
% cd ..
% cd 5Xi
% save PT.txt -ascii Ept_Xi;
% fid = fopen('jobNum.txt','wt'); 
% fprintf(fid,'%g\n',length(Ept_Xi));       
% fclose(fid);
% cd ..