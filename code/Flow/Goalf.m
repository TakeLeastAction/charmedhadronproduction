function y=Goalf(param)
global   Ept_phi Espt_phi Eerr_phi  Ept_Omega  Espt_Omega  Eerr_Omega  Ept_p  Espt_p  Eerr_p  Ept_Lambda  Espt_Lambda Eerr_Lambda  Ept_Xi  Espt_Xi Eerr_Xi
global     IT
clc


 load ExpData  
 [Ept_phi, Espt_phi,Eerr_phi, Ept_Omega, Espt_Omega,Eerr_Omega,Ept_p, Espt_p,Eerr_p, Ept_Lambda, Espt_Lambda,Eerr_Lambda, Ept_Xi, Espt_Xi,Eerr_Xi]=Initial();


[SPT_phi,SPT_Omega,SPT_p,SPT_Lambda,SPT_Xi]=Func1(param);
SPT_Xi = SPT_Xi*3;%% Important!
SPT_p = SPT_p*5;
SPT_Lambda = SPT_Lambda*7.44;
k=0;
while k==0
    try
    y1=mean((Espt_phi-SPT_phi).^2./Eerr_phi.^2);  
    y2=mean((Espt_Omega-SPT_Omega).^2./Eerr_Omega.^2);  
    y3=mean((Espt_p-SPT_p).^2./Eerr_p.^2);  
    y4=mean((Espt_Lambda-SPT_Lambda).^2./Eerr_Lambda.^2);  
    y5=mean((Espt_Xi-SPT_Xi).^2./Eerr_Xi.^2);  
    k=1;
    catch
        load ExpData;
    end
end

%y=(y2+y4+y5)/3;
y=mean([y1 y2 y3 y4 y5]);

%hold on;plot()

fid=fopen('process.txt','a');
fprintf(fid,'%d\n',[]);
fprintf(fid,'Iteration number: ');
fprintf(fid,'%d\n',IT);
fprintf(fid,'parameters : ');
fprintf(fid,'%d\n',[]);
fprintf(fid,'%d\n', param);
fprintf(fid','Error of phi is :');
fprintf(fid,'%d',y1);
fprintf(fid,'%d\n',[]);
fprintf(fid','Error of Omega is :');
fprintf(fid,'%d',y2);
fprintf(fid,'%d\n',[]);
fprintf(fid','Error of proton is :');
fprintf(fid,'%d',y3);
fprintf(fid,'%d\n',[]);
fprintf(fid','Error of Lambda is :');
fprintf(fid,'%d',y4);
fprintf(fid,'%d\n',[]);
fprintf(fid','Error of xi is :');
fprintf(fid,'%d',y5);
fprintf(fid,'%d\n',[]);
fprintf(fid','Error of objFun is :');
fprintf(fid,'%d',y);

fclose(fid);

disp( 'Iteration number:')
disp( IT)

% subject1=strcat('It=',num2str(IT),';objF is',num2str(y),';','P=',...
%     num2str(param.'),';','phi=',num2str(y1),';Om=',num2str(y2),';p=',...
%     num2str(y3));
% subject2=strcat('It=',num2str(IT),';Lam=',num2str(y4),';xi=',num2str(y5));
% content='Kiss while your lips are red';
% 
% s=toc;
% if s>1800
%     mail2me(subject1,content);
%     mail2me(subject2,content);
%     tic
% end

subject=strcat('Hi Master: This is the fitting of spectrum of LHC. The iteration number is ',num2str(IT),', objF is',num2str(y),', param is',...
    num2str(param.'),', phi is',num2str(y1),', Omega is',num2str(y2),', proton is',...
    num2str(y3),', Lamda is',num2str(y4),', xi is',num2str(y5));

content='   Kiss while your lips are red!';
text=strcat(subject,content);

s=toc;
if s>600
    %mail2me(subject,content);
    tic
end

IT=IT+1;


