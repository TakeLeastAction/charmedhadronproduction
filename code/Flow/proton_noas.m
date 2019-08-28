function SPT=proton_noas(qm0,PT1,T,beta0,tao0,R0,as,c1,c2,s2,c_v1,c_v3,c_v4,s4,B)

%global c_v1 c_v4 c_v3 s4

% if T<=0.165
% tao0=exp(exp(tao0));
% end
tao0=exp(tao0);
NN=100;
hba=0.19733;
%PT=3;
%qm0=0.938;
%R0=16.0;
%tao0=113.2;  %/6.0*(16/16)^2
%y=0;
r_max=1;
r_min=0;
Fais_max=2*pi;
Fais_min=0;
%Rz_max=tao0*sinh(1.5);
%Rz_min=-Rz_max;
Rba=linspace(r_min,r_max,NN);
Fais=linspace(Fais_min,Fais_max,NN+1);
%Rz1=linspace(Rz_min,Rz_max,NN);
%tao1=linspace(tao1min,tao1max,NN);
[rba,fais]=meshgrid(Rba,Fais);     
for k=1:length(PT1)
    PT=PT1(k);
    qmT=sqrt(qm0^2+PT^2);
    %Pz=qmT*sinh(y);
      %s2=0.0 ;
      %T=0.1107 ;
      %c1=0.0 ;
      %c2=3.0 ;
      %as=0.01;
      %beta0=0.9801 ;      
      %N=1;
      %B=1 ;
     % FB=0 ;
      %tf1=0.0;
      %tf2=0.0;
      %NL=0   ;     
      rbound=R0*(1+s2*cos(2*fais)+s4*cos(4*fais));
      rbound_1=R0*(-2*s2*sin(2*fais)-4*s4*sin(4*fais));
      tgfb=-(-rbound.*tan(fais)+rbound_1)./(rbound+tan(fais).*rbound_1);
      
      
      faib=atan(tgfb);%
	  p=fais>pi/2; 
      faib(p)=faib(p)+pi;
      p=fais>3*pi/2; 
      faib(p)=faib(p)+pi;
      %faib=fais;
      
     %[fais,r] =cart2pol(Rx,Ry);  
     %[faip,PT] =cart2pol(Px,Py);
      int1=zeros(NN+1,NN);
      %y=asinh(Pz/qmT);   
      if as==0
          wmga=1;
      else
      wmga=1./(1+exp((rba-1)/(as+eps)));
      end
      %ro=beta0*rba.*(1+c1*exp(-PT/c2)*cos(2*faib));
      vT = PT./qm0;
      r=R0*(1+s2)*(1-s2)*rba./sqrt(1+s2^2+2*s2*((sin(fais)).^2-(cos(fais)).^2));%pujie 
      ro=beta0*rba.*(1+exp(-r.^2.0./c2^2).*(c1*cos(2*faib)+c_v1*cos(faib)+c_v3*cos(3*faib)+c_v4*cos(4*faib)));%pujie
      %20181011pujie%ro=beta0*rba.*(1+exp(-c2*rba.^2).*(c1*cos(2*faib)+c_v1*cos(faib)+c_v3*cos(3*faib)+c_v4*cos(4*faib)));
      %ro=beta0*rba.*(1+exp(-vT/c2)*(c1*cos(2*faib)+c_v1*cos(faib)+c_v3*cos(3*faib)+c_v4*cos(4*faib)));
      afa=PT/T.*sinh(ro);
      beta=qmT/T.*cosh(ro);

      Btheta=1./rbound; 
      
      
      if B==0
      int1(:,:)=2*(2*pi*0.19733)^(-3).*tao0*qmT*2*besselk(1,beta).*wmga.*rba.*besseli(0,afa)./Btheta.^2;
      else if B==1
              int1(:,:)=2*(2*pi*0.19733)^(-3).*tao0*qmT*2*besselk(1,beta).*wmga.*rba.*besseli(1,afa)./Btheta.^2.*cos(faib); 
      else if B==2
             int1(:,:)=2*(2*pi*0.19733)^(-3).*tao0*qmT*2*besselk(1,beta).*wmga.*rba.*besseli(2,afa)./Btheta.^2.*cos(2*faib); 
          else if B==3
                  int1(:,:)=2*(2*pi*0.19733)^(-3).*tao0*qmT*2*besselk(1,beta).*wmga.*rba.*besseli(3,afa)./Btheta.^2.*cos(3*faib);
              else if B==4
                      int1(:,:)=2*(2*pi*0.19733)^(-3).*tao0*qmT*2*besselk(1,beta).*wmga.*rba.*besseli(4,afa)./Btheta.^2.*cos(4*faib); 
                  else if B>4
                        int1(:,:)=2*(2*pi*0.19733)^(-3).*tao0*qmT*2*besselk(1,beta).*wmga.*rba.*besseli(B,afa)./Btheta.^2.*cos(B*faib);   
                      end
                  end
              end
          end
          end
      end
      %int2=zeros(NN);
      int2=trapz(Fais,int1);
      SPT(k)=trapz(Rba,int2);
      %figure;plot(Rba,int2,'-r.')
end

      %figure;plot(R,int2,'-r.');xlabel('r');ylabel('f')