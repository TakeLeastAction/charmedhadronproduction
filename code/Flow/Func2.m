function SPT=Func2(qm0,PT1,T,beta0,xi)
% dN/2*pi*ptdptd(yita)   proton_noas
% xi is the fugacity

s2=0.0;c1=0.0;c2=0.0;faip=0.0;as=0.0;
R0=15.6;
tao0=15.2;

NN=100;
hbac=0.19733;
%faip=0;
R_max=R0;
R_min=0;
Fais_max=2*pi;
Fais_min=0;
R=linspace(R_min,R_max,NN);
Fais=linspace(Fais_min,Fais_max,NN);
[r,fais]=meshgrid(R,Fais);

for k=1:length(PT1)
    PT=PT1(k);
    Px=PT*cos(faip);
    Py=PT*sin(faip);
    qmT=sqrt(qm0^2+PT^2);
    Rx0=R0*(1+s2);
    Ry0=R0*(1-s2); 
     %[fais,r] =cart2pol(Rx,Ry);  
     %[faip,PT] =cart2pol(Px,Py);
    faib=fais;
    int1=zeros(NN,NN);
      %y=asinh(Pz/qmT);   
    rba=sqrt((r.*cos(fais)./Rx0).^2+(r.*sin(fais)./Ry0).^2);
    wmga=1./(1+exp((rba-1)/(as+eps)));
    ro=beta0*rba.*(1+c1*exp(-PT/c2)*cos(2*faib));
    afa=PT/T.*sinh(ro);
    beta=qmT/T.*cosh(ro);
    p1=afa.*cos(faib-faip);
    f=exp(p1);
    int1(:,:)=tao0*qmT*2*besselk(1,beta).*wmga.*r.*f.*2*(2*pi*hbac)^(-3) ;       
    int2=trapz(Fais,int1);
    SPT(k)=trapz(R,int2);
end
SPT=SPT.'*xi;