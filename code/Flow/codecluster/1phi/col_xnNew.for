clw===========================================================
      external coles
      dimension param(50)
      dimension pt_phi(100)
      common/input/T,R0,tao0,taod,yidu,beta0,c1,c2,s2,ell_1,
     c ell_3,ell_4,s4,qmu,qms,qmc,yidu_u,yidu_anti_u,yidu_s,
     c yidu_anti_s, yidu_c,yidu_anti_c,freq
      open(5,file='input.txt')
      read(5,*)T
      read(5,*)beta0
      read(5,*)R0
      read(5,*)tao0
      read(5,*)taod
      read(5,*)qmu
      read(5,*)qms
      read(5,*)qmc
      read(5,*)freq
      read(5,*)yidu_u
      read(5,*)yidu_anti_u
      read(5,*)yidu_s     
      read(5,*)yidu_anti_s
      read(5,*)yidu_c
      read(5,*)yidu_anti_c   
      read(5,*)c1
      read(5,*)c2
      read(5,*)s2
      read(5,*)ell_1
      read(5,*)ell_3
      read(5,*)ell_4
      read(5,*)s4
      open(6,file='circle.txt')
      read(6,*)circle
      open(7,file='jobNum.txt')
      read(7,*)jobNum
      open(8,file='PT.txt')
      do 110 ij=1,jobNum
             read(8,*)pt_phi(ij)
110   continue

      do 02 j=circle,circle   ! 横动量
      do 03 k=1,1    ! 角度 101 central 
      do 01 i=1,50    ! 初始化
      param(i)=0.0
 01   continue
      param(1)=pt_phi(j)!(j-1)*0.875/4.0  !(j-1)/1.     !           ! 横动量
      param(2)=3.1415926/50.0*(k-1) ! 角度
      param(3)=2            ! 粒子数
      param(4)=qms !1.115 !
      param(5)=qms !1.115 !
      param(6)=param(5)
      param(7)=param(5)
      param(8)=param(5)
      param(9)=param(5)
      n=param(3)
      call vegas(coles,0.0001,n*7-3,500000,100,0,0,param,yout,
     c yout_err ,yout_sd)
      open(40,file="phi.txt")  !Lam_H5_pm2_10_15
      write(40,*)yout
      open(41,file="phi_err.txt")
      !Lam_H5_pm2_10_15
      write(41,*)yout_err
      open(42,file="phi_sd.txt") 
       !Lam_H5_pm2_10_15
      write(42,*)yout_sd
 03   continue
 02   continue
      
	
	stop
	end

clw----------------------1/(sqrt(2*pi)*Pc*Rc)**3-------------------------------------
	function coles(x,param)

	dimension x(40) , param(10),nn(40)
      common/input/T,R0,tao0,taod,yidu,beta0,c1,c2,s2,ell_1,
     c ell_3,ell_4,s4,qmu,qms,qmc,yidu_u,yidu_anti_u,yidu_s,
     c yidu_anti_s, yidu_c,yidu_anti_c,freq
	data pi/3.1415926/
	data RMS/0.59/ 
	data deta1/3.0/
      data deta2/3.0/
      data deta3/3.0/
      !data taod/4.0/
	data hba/0.19733/
      dimension fun(param(3)), qm0(param(3))
      dimension qmt(param(3)),qmiu(param(3)),qmt_1(param(3))
      dimension Pb(param(3)-1),Pc(param(3)-1)
      dimension tao(param(3)),r(param(3),3),pop(param(3),3)
      dimension rel_r(param(3),3),rel_p(param(3),3)
      dimension qJaco(param(3),param(3)),tri(param(3),param(3))
      dimension qJ_1(param(3),param(3)),qJ_t(param(3),param(3))
      dimension qJ_1t(param(3),param(3)),ftao(param(3))
      dimension yf(param(3)),rel_rnew(param(3),3)
      dimension rel_pnew(param(3),3),rel_rnew1(param(3),3)
      dimension rel_rmax(param(3),3),rel_pmax(param(3),3)
      dimension rel_rmin(param(3),3),rel_pmin(param(3),3)
      dimension tao_max(param(3)),tao_min(param(3))
      dimension qr(param(3),3),qp(param(3),3),qLr(param(3),3)
      dimension r_cy(param(3),3),p_cy(param(3),3),qLp(param(3),3)
      dimension rLE(param(3)),rLx(param(3)),rLy(param(3)),rLz(param(3))
      dimension ppE(param(3)),ppx(param(3)),ppy(param(3)),ppz(param(3))
      dimension yta(param(3)),y_y(param(3)),qmt1(param(3)),q_m(param(3))
      dimension faib(param(3)),ro(param(3)),rba(param(3))
      dimension afa(param(3)),beta(param(3)),wmga(param(3))
      dimension sig2(param(3)-1)

      N=param(3)
      col_fact=8.0 
      faip=param(2)    
      PT=param(1)
      qmtotoal=0.0

      do 801 i=1,N
      qm0(i)=param(i+3)
801   continue
      qmt(1)=qm0(1)
      qmt_1(1)=1./qm0(1)
      do 901 i=1,N-1
      qmt(i+1)=qmt(i)+qm0(i+1)
      qmt_1(i+1)=qmt_1(i)+1./qm0(i+1)
901   continue


      qmtotoal=qmt(N)

      do 802 i=1,N-1
      qmiu(i)=(i+1.)/(i+0.)*(1.0/qmt(i)+1.0/qm0(i+1))**(-1.)
802   continue

      !wmga_i=2.0*(N+0.)/3.*RMS**2.0*1.0/(qmt_1(N)-(N+0.)/qmtotoal)
      wmga_i = freq ! frequency of hamornic oscillator
      do 803 i=1,N-1
      sig2(i)=wmga_i/qmiu(i)
      Pb(i)=sqrt(sig2(i))
      Pc(i)=hba/Pb(i)
803   continue

      y=0.0
      !!R0=16.0 !param(1) !5;
      !s2=-0.0818
      !yidu=1.0 !6.84 !7.597  !7.913   ! 
      !tao0=3.446/yidu*(16.0/R0)**2.0
      rel_rmax(1,1)=R0*(1.0+abs(s2)+abs(s4))
      rel_rmin(1,1)=-rel_rmax(1,1)
      R_x=(rel_rmax(1,1)-rel_rmin(1,1))*x(1)+rel_rmin(1,1) 
      rel_rmax(1,2)=sqrt(rel_rmax(1,1)**2-R_x**2)
      rel_rmin(1,2)=-rel_rmax(1,2)
      rel_rmax(1,3)=tao0*sinh(1.1)  
      rel_rmin(1,3)=-rel_rmax(1,3)
      rel_r(1,1)=R_x
      rel_r(1,2)=(rel_rmax(1,2)-rel_rmin(1,2))*x(2)+rel_rmin(1,2) 
      rel_r(1,3)=(rel_rmax(1,3)-rel_rmin(1,3))*x(3)+rel_rmin(1,3) 

      do 804 i=1,N-1
      rel_rmax(i+1,1)=1/2.000000+atan(deta1*Pb(i))/pi
      rel_rmin(i+1,1)=1/2.000000-atan(deta1*Pb(i))/pi
      rel_rmax(i+1,2)=rel_rmax(i+1,1)
      rel_rmax(i+1,3)=rel_rmax(i+1,1)
      rel_rmin(i+1,2)=rel_rmin(i+1,1)
      rel_rmin(i+1,3)=rel_rmin(i+1,1)
      rel_r(i+1,1)=(rel_rmax(i+1,1)-rel_rmin(i+1,1))
     c *x(i*3+1)+rel_rmin(i+1,1) 
      rel_r(i+1,2)=(rel_rmax(i+1,2)-rel_rmin(i+1,2))
     c *x(i*3+2)+rel_rmin(i+1,2) 
      rel_r(i+1,3)=(rel_rmax(i+1,3)-rel_rmin(i+1,3))
     c *x(i*3+3)+rel_rmin(i+1,3) 


       rel_r(i+1,1)=tan(pi*(rel_r(i+1,1)-1./2.))
       rel_r(i+1,2)=tan(pi*(rel_r(i+1,2)-1./2.))
       rel_r(i+1,3)=tan(pi*(rel_r(i+1,3)-1./2.))
804   continue


  
      do 805 i=1,N-1
      rel_pmax(i+1,1)=1/2.000000+atan(deta2*Pc(i))/pi
      rel_pmin(i+1,1)=1/2.000000-atan(deta2*Pc(i))/pi
      rel_pmax(i+1,2)=rel_pmax(i+1,1)
      rel_pmax(i+1,3)=rel_pmax(i+1,1)
      rel_pmin(i+1,2)=rel_pmin(i+1,1)
      rel_pmin(i+1,3)=rel_pmin(i+1,1)
      rel_p(i+1,1)=(rel_pmax(i+1,1)-rel_pmin(i+1,1))
     c *x(i*3-2+N*3)+rel_pmin(i+1,1) 
      rel_p(i+1,2)=(rel_pmax(i+1,2)-rel_pmin(i+1,2))
     c *x(i*3-1+N*3)+rel_pmin(i+1,2) 
      rel_p(i+1,3)=(rel_pmax(i+1,3)-rel_pmin(i+1,3))
     c *x(i*3+N*3)+rel_pmin(i+1,3) 

       rel_p(i+1,1)=tan(pi*(rel_p(i+1,1)-1./2.))
       rel_p(i+1,2)=tan(pi*(rel_p(i+1,2)-1./2.))
       rel_p(i+1,3)=tan(pi*(rel_p(i+1,3)-1./2.))
805   continue

      do 806 i=1,N
      tao_max(i)=tao0+deta3*taod
      tao_min(i)=tao0-deta3*taod
806   continue
      do 906 i=1,N
      tao(i)=(tao_max(i)-tao_min(i))*x(6*N-3+i)+tao_min(i)
906   continue      

      do 807 i=1,N
      qJaco(1,i)=qm0(i)/qmt(N)
807   continue
 
      do 808 i=1,N
      do 809 j=1,N
      tri(i,j)=0.
809   continue
808   continue

      do 810 i=1,N-1
      tri(i,i+1)=1.
810   continue
   

      do 812 i=1,N-1
      do 811 j=1,N
      if (i.ge.j) then
       temp=qm0(j)
      else 
      temp=0.
      end if
      qJaco(i+1,j)=sqrt((i+0.0)/(i+1.0))*(temp/qmt(i)-tri(i,j))
811   continue
812   continue

      do 813 i=1,N
      qJ_1(i,1)=1.
813   continue

      do 814 i=1,N
      do 815 j=2,N
      if (i.gt.j) then
       qJ_1(i,j)=0.
      else if(i.eq.j)then
           qJ_1(i,j)= -sqrt((j+0.)/(j-1.))*(qmt(j-1)/qmt(j))
           else
           qJ_1(i,j)= sqrt((j+0.)/(j-1.))*(qm0(j)/qmt(j))
      end if
815   continue
814   continue


      do 816 i=1,N
      do 817 j=1,N
      qJ_t(i,j)=qJaco(j,i)
      qJ_1t(i,j)=qJ_1(j,i)
817   continue
816   continue


      P_x=PT*cos(faip)
      P_y=PT*sin(faip)
      qmTot=sqrt(qmt(N)**2+PT**2)
      P_z=qmTot*sinh(y)

      rel_p(1,1)=P_x
      rel_p(1,2)=P_y
      rel_p(1,3)=P_z


      do 818 i=1,N
      do 819 j=1,3
      qr(i,j)=0.
      qp(i,j)=0.
      do 820 k=1,N
      qr(i,j)=qr(i,j)+qJ_1(i,k)*rel_r(k,j)
      qp(i,j)=qp(i,j)+qJ_t(i,k)*rel_p(k,j)
820   continue
819   continue
818   continue

      do 821 i=1,N
      call Cychange(qr(i,1),qr(i,2),r_cy(i,1),r_cy(i,2))
      call Cychange(qp(i,1),qp(i,2),p_cy(i,1),p_cy(i,2))
821   continue
      N_yta=0
      do 822 i=1,N
      qmt1(i)=sqrt(qm0(i)**2.+p_cy(i,1)**2.)
      q_m(i)=sqrt(qm0(i)**2.+p_cy(i,1)**2.+qp(i,3)**2.)
      yta(i)=asinh(qr(i,3)/tao(i))
     
      if (abs(tao(i)).lt.1e-10)then
       N_yta=1
       write (*,*) i,tao(i)
      end if
      y_y(i)=asinh(qp(i,3)/qmt1(i))
822   continue     
      


      !T=0.1115 !175 !param(3)!0.175;              %   GeV
      !c1=0.1111 !0.43 !param(4)!0.43;            %
      !c2=4.0188 !param1(1)!0.85;          %
      !as=0.00043 !param1(2)!0;             %
      !beta0=0.9185 !0.55 !param1(4)!0.55;
      !q=1.0095

      B=1 !param2(3)!1;  
      FB=0 !param2(4)!1;
      tf1=0.0
      tf2=0.0
      NL=0
      !Rx=R0*(1+s2)
      !Ry=R0*(1-s2)
      Bf=1
      i=1
      if (B.eq.1)then
        Bf=1
        else if(B.eq.2)then
          Bf=cos(2*faip)
       end if

      do 823 i=1,N
      rbound=R0*(1+s2*cos(2*r_cy(i,2))+s4*cos(4*r_cy(i,2)))   
      ! fais=r_cy(i,2)
      rbound1=R0*(-2*s2*sin(2*r_cy(i,2))-4*s4*sin(4*r_cy(i,2)))
      tgfb=-(-rbound*tan(r_cy(i,2))+rbound1)
     c /(rbound+tan(r_cy(i,2))*rbound1)
      faib(i)=atan(tgfb)
      if(r_cy(i,2)>pi/2)then
         faib(i)=faib(i)+pi
      end if
      if(r_cy(i,2)>3*pi/2)then
      faib(i)=faib(i)+pi
      end if
823   continue         
      do 824 i=1,N
      rba(i)=sqrt(qr(i,1)**2.+qr(i,2)**2.)
     c /(R0*(1+s2*cos(2*r_cy(i,2))+s4*cos(4*r_cy(i,2)))) 
       if (rba(i)>1.0)then
                    wmga(i)=0.0
        else if(rba(i)<1.0)then
                    wmga(i)=1.0
       end if
      !wmga1=1/(1+exp((rba1-1)/(as+0.0000000001)))
      ro(i)=beta0*rba(i)*(1+exp(-p_cy(i,1)/c2)*(c1*cos(2*faib(i))
     c +ell_1*cos(faib(i))+ell_3*cos(3*faib(i))+ell_4*cos(4*faib(i))))
      afa(i)=p_cy(i,1)/T*sinh(ro(i))
      beta(i)=qmt1(i)/T*cosh(ro(i))
824    continue


        do 8 i=1,N
        fun(i)=exp(afa(i)*cos(faib(i)-p_cy(i,2))-beta(i)*
     c  (cosh(yta(i)-y_y(i))))
 8       continue
        do 825 i=1,N
        ftao(i)=1/taod/sqrt(2*pi)*exp(-(tao(i)-tao0)**2/2.0/taod**2.0)
        yf(i)=(cosh(y_y(i)-yta(i)))/cosh(y_y(i))/cosh(yta(i))
     c *wmga(i)*fun(i)*2*(2*pi*0.19733)**(-3)*ftao(i) !10.0**3*
 825       continue

      qmtotal=0.
        do 826 i=1,N
        qmtotal=qmtotal+q_m(i)
 826  continue

      bex=P_x/qmtotal
      bey=P_y/qmtotal
      bez=P_z/qmtotal
      !b0=sqrt(1+bex**2+bey**2+bez**2)
      !bex=bex/b0
      !bey=bey/b0
      !bez=bez/b0
      
      !before transformed     
      do 827 i=1,N
      CALL lorentz(tao(i)*cosh(yta(i)),qr(i,1),qr(i,2),
     c qr(i,3),
     c bex,bey,bez,
     c rLE(i),rLx(i),rLy(i),rLz(i))
      CALL lorentz(qmt1(i)*cosh(y_y(i)),qp(i,1),qp(i,2),
     cqp(i,3),
     cbex,bey,bez,
     c ppE(i),ppx(i),ppy(i),ppz(i))
 827  continue


      do  914 i=1,N
      do  915 j=1,3
      if (j.eq.1) then
      qLr(i,j)=rLx(i)
      qLp(i,j)=ppx(i)
      else if (j.eq.2) then 
      qLr(i,j)=rLy(i)
      qLp(i,j)=ppy(i)
      else 
      qLr(i,j)=rLz(i)
      qLp(i,j)=ppz(i)
      end if
915   continue
914   continue
     



      do 828 i=1,N
      do 829 j=1,3
      rel_rnew(i,j)=0.
      rel_pnew(i,j)=0.
      do 830 k=1,N
      rel_rnew(i,j)=rel_rnew(i,j)+qJaco(i,k)*qLr(k,j)
      rel_pnew(i,j)=rel_pnew(i,j)+qJ_1t(i,k)*qLp(k,j)
830   continue
829   continue
828   continue     
     
     
     
     
     
     
      

      do 910 i=1,N-1
      if (i.eq.1) then
      s=max(rLE(i),rLE(i+1))
      else
      s=max(s,rLE(i+1))
      end if
910   continue
      TTT=s
      do 911 i=1,N
       rLx(i)=rLx(i) +(TTT-rLE(i))*ppx(i)/ppE(i)
       rLy(i)=rLy(i) +(TTT-rLE(i))*ppy(i)/ppE(i)
       rLz(i)=rLz(i) +(TTT-rLE(i))*ppz(i)/ppE(i)
911   continue
 


      do  912 i=1,N
      do  913 j=1,3
      if (j.eq.1) then
      qLr(i,j)=rLx(i)
      else if (j.eq.2) then 
      qLr(i,j)=rLy(i)
      else 
      qLr(i,j)=rLz(i)
      end if

913   continue
912   continue


      
      
      do 831 i=1,N
      do 832 j=1,3
      rel_rnew1(i,j)=0.
      do 833 k=1,N
      rel_rnew1(i,j)=rel_rnew1(i,j)+qJaco(i,k)*qLr(k,j)
833   continue
832   continue
831   continue    

      

     
      
        wigner_q=0.
        wigner_p=0.
        do 834 i=1,N-1
        wigner_q=wigner_q+(rel_r(i+1,1)**2.+rel_r(i+1,2)**2.
     c +rel_r(i+1,3)**2.)/Pb(i)**2.
        wigner_p=wigner_p+(rel_p(i+1,1)**2.+rel_p(i+1,2)**2.
     c +rel_p(i+1,3)**2.)/Pc(i)**2.
 834    continue

         wigner=wigner_q+wigner_p
        fw1=col_fact**(N-1)*exp(-wigner)

        wigner_qnew=0.
        wigner_pnew=0.
        do 837 i=1,N-1
        wigner_qnew=wigner_qnew+(rel_rnew1(i+1,1)**2.+rel_rnew1(i+1,2)
     c **2.+rel_rnew1(i+1,3)**2.)/Pb(i)**2.
        wigner_pnew=wigner_pnew+(rel_pnew(i+1,1)**2.+rel_pnew(i+1,2)**2.
     c +rel_pnew(i+1,3)**2.)/Pc(i)**2.
 837    continue

        wigner=wigner_qnew+wigner_pnew
        fw2=col_fact**(N-1)*exp(-wigner)


      if ((N.eq.2).and.(qm0(1).eq.0.938)) then
      X1=rLx(1)
      X2=rLx(2)
      Y1=rLy(1)
      Y2=rLy(2)
      Z1=rLz(1)
      Z2=rLz(2)   

      PX1=ppx(1)
      PX2=ppx(2)
      PY1=ppy(1)
      PY2=ppy(2)
      PZ1=ppz(1)
      PZ2=ppz(2)

      DR =SQRT( (X1-X2)**2 + (Y1-Y2)**2 + (Z1-Z2)**2 )
      DP =0.5*SQRT( (PX1-PX2)**2 + (PY1-PY2)**2 + (PZ1-PZ2)**2 )
      rdotk=((X1-X2)*(PX1-PX2)+(Y1-Y2)*(PY1-PY2)+(Z1-Z2)*(PZ1-PZ2))
     &/hba/2.

      fw3 = HulthenWG15(DR,DP,rdotk)      
      end if 

        
        acu_tao=1.
        acu_f=1.
        do 835 i=1,N
        acu_tao=(tao_max(i)-tao_min(i))*acu_tao
        acu_f=acu_f*yf(i)
 835    continue

         acu_q=1.
         acu_p=1.
        do 836 i=1,N-1
        acu_q=(1+rel_r(i+1,1)**2.)*(1+rel_r(i+1,2)**2.)
     c  *(1+rel_r(i+1,3)**2.)*(2.000*atan(deta1*Pb(i)))**3.*acu_q
        acu_p=(1+rel_p(i+1,1)**2.)*(1+rel_p(i+1,2)**2.)
     c  *(1+rel_p(i+1,3)**2.)*(2.000*atan(deta2*Pc(i)))**3.*acu_p
 836    continue




      S=(rel_rmax(1,1)-rel_rmin(1,1))*(rel_rmax(1,2)-rel_rmin(1,2))
     c *(rel_rmax(1,3)-rel_rmin(1,3))*acu_q*acu_p*acu_tao

   
     
      coles=qmTot*cosh(y)*acu_f*Bf*S*fw2*yidu_s*yidu_anti_s**N 
      s=yidu**N
      if (N_yta.eq.1)then
      coles=0.0
      end if   
      !if (isnan(coles))then
      !coles=0.0
      !end if 
      !write(*,*) coles

	return
	end
	


CLW\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\C

       function HulthenWG15(r,q,rdotk)

CLW\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\C
      dimension A(15),C(15)
      data A/3.496647458963869E-001,2.514620979752955E-004,
     & 1.062940427959510E-001,3.229472385579899E-002,
     & 4.890471441137741E-002,1.722798854593734E-001,
     & 7.742048546701388E-002,1.726701557237221E-002,
     & 1.482683976833659E-004,4.198161133325294E-002,
     & 1.854191668060074E-001,7.355493587201123E-003,
     & 4.621519814465522E-002,1.158256453709526E-002,
     & 1.494576593001586E-001/
      data C/1.579571884334015E-002,1.988218908245876,
     & 9.276209379806834E-001,2.592426417692896,
     & 2.899018075383141E-001,8.997928431726522E-002,
     & 1.933532774200236E-001,4.896043560931306E-001,
     & 1.998110600951886E-001,4.707391003810082E-001,
     & 3.942925929143137E-002,2.179211335552845E-001,
     & 9.759425154944455E-002,14.463899746746620,
     & 1.801166409712259E-001/

      N=15
      q=q/0.19733
c	rdotk=0.

      T1=0.
	DO I=1,N
	T1=T1+A(I)*A(I)*EXP(-2.*C(I)*r*r-0.5*q*q/C(I))
	END DO

      T2=0.
	DO I=1,N-1
	   DO J=I+1,N
	T2=T2+A(I)*A(J)*(4.*ABS(C(I)*C(J))/(C(I)+C(J))**2)**0.75
     &*EXP(-4.*C(I)*C(J)/(C(I)+C(J))*r*r-q*q/(C(I)+C(J)))
     &*COS(2.*(C(I)-C(J))/(C(I)+C(J))*rdotk)
	   END DO
	END DO
     
      HulthenWG15=8.*T1+16*T2
	
      RETURN
	END     	
	
	subroutine Cychange(rx,ry,r,fai)
	r=sqrt(rx**2+ry**2)
	fai=acos(rx/r)
	if (ry<0) then
	fai=2*3.1415926-fai
	end if
	return 
	end
	
	
      subroutine lorentz(energy, px, py, pz, bex, bey, bez,
     cenergy_new, px_new, py_new, pz_new)
CLW\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\C
CLW: Frame F    : energy, px, py, pz 
CLW: Frame F_new: energy_new, px_new, py_new, pz_new
CLW:           Frame F_new moves with velocity bex, bey, bez 
CLW: with respect to Frame F. 
CLW:              (energy, px, py, pz<===>time,x,y,z)
c     add in a cut for beta2 to prevent gam to be nan (infinity)

      !implicit real*8 (a-h, o-z)

CLW      common /lor/ energy_new, px_new, py_new, pz_new

      beta2 = bex ** 2 + bey ** 2 + bez ** 2
      if (beta2 .eq. 0d0) then
         energy_new = energy
         px_new = px
         py_new = py
         pz_new = pz
      else
         if (beta2 .gt. 0.999999999999999d0) then
            beta2 = 0.999999999999999d0
            print *,'beta2=0.999999999999999'
         end if
         b0 = dsqrt(1.d0+beta2)
         
         energy_new = (b0 * energy - bex * px - bey * py - bez * pz)
         px_new = -  bex * energy + (1.d0 
     c        + 1.d0/(1+b0) * bex ** 2 ) * px
     c        + 1.d0/(1+b0)* bex * bey * py
     c        + 1.d0/(1+b0)* bex * bez * pz  
       !c1= b0**2/(1+b0)
       !c2=(b0-1)/beta2
       
         py_new = - bey * energy 
     c        + 1.d0/(1+b0) * bex * bey * px
     c        + (1.d0 + 1.d0/(1+b0) * bey ** 2 ) * py
     c        + 1.d0/(1+b0) * bey * bez * pz         
         pz_new = - bez * energy
     c        +  1.d0/(1+b0)* bex * bez  * px
     c        +1.d0/(1+b0)* bey * bez * py
     &        + (1.d0 + 1.d0/(1+b0) * bez ** 2) * pz    
      endif

      return
      end
	

c====================================================================
c   01/12/84 602042308	member name  vegas    (s)	    fortran
c ===================================================================
c
c     vegas version including hbook calls possibility,
c     simplified but more clear and equally efficient.
c
c			       xd.zhang imp-92
c ===================================================================
c
c     fxn   = function to be integrated/mapped
c     acc   = relative accuracy requested
c     ndim  = # dimensions
c     ncall = maximum total # of calls to the function per iteration
c     itmx  = maximum # of iterations allowed
c     nprn  = printout level:
c	      =0  only final results
c	     >=1  additionnally inf. about input parameters
c	     >=2  additionnally inf. about accumulated values
c		   per iteration.
c	     >=3  additionnally inf. about partial values
c		   per iteration.
c	     >=5  additionnally inf. about final bin distribution
c		   (numerical mapping).
c     igraph= histogrammation level:
c	      =0  no histograms at all
c	      =1  only statistics about the integration
c	      =10 only histograms defined but not statistics
c
c
c ===================================================================
      subroutine vegas(fxn,acc,ndim,ncall,itmx,nprn,igraph,param,yout,
     c yout_err,yout_sd)
c     implicit real*8(a-h,o-z)
      common/result/it,chi2a,avgi,sd,ti,tsi 
      common/gas/deint
      common/inparm/itmx0,alph
      common/pp/wl,r1,iseed,np
c     common/ee/iseed,ird
      common/outmap/nd,ndim0,xi(50,50)
      dimension x(50),xin(50),r(50),ia(50),d(50,50),param(12)
      data alph0/1.5/,init/0/
      iseed=999999999
      nd=50
c===========================================================
c a))  initializing some variables
c===========================================================
      itmx0=itmx
      ndim0=ndim
      if(alph.eq.0.)alph=alph0
      calls=ncall
      xnd=nd
      ndm=nd-1
c     if(igraph.ne.0)call inbook(0,fun,weight,igraph)
c.............................................
c  initializing cummulative variables
c.............................................
      it=0
      si=0.
      si2=0.
      swgt=0.
      schi=0.
      scalls=0.
      nzerot=0
      fmaxt=0.
      fmint=1.e30
c.............................................
c  defining the initial intervals distribution
c.............................................
      rc=1./xnd
      do 7 j=1,ndim
      xi(nd,j)=1.
      dr=0.
      do 7 i=1,ndm
      dr=dr+rc
      xi(i,j)=dr
7     continue
      if(init.eq.0)print 407
      init=1
      if(nprn.ge.1)print 290,ndim,ncall,itmx,nd
c===========================================================
c b))  iterations loop
c===========================================================
9     it=it+1
c.............................................
c  initializing iteration variables
c.............................................
      ti=0.
      sfun2=0.
      nzero=0
      fmax=0.
      fmin=1.e30
      do 10 j=1,ndim
      do 10 i=1,nd
      d(i,j)=0.
10    continue
c     if(igraph.ne.0)call rebook(0,fun,weight,igraph)
      do 11 jj=1,ncall
      wgt=1.
c.............................................
c  computing the point position
c.............................................
      do 15 j=1,ndim
      xn=ran(iseed)*xnd+1.
      ia(j)=xn
      xim1=0.
      if(ia(j).gt.1)xim1=xi(ia(j)-1,j)
      xo=xi(ia(j),j)-xim1
      x(j)=xim1+(xn-ia(j))*xo
      if (x(j).le.1e-20)x(j)=1e-20
      wgt=wgt*xo*xnd
15    continue
c.............................................
c  computing the function value
c.............................................

      fun=fxn(x,param)
      !if (it.eq.18) open(30,file="re.txt")
      !if (it.eq.18) then
          !if (jj.eq.21769) write(30,*)x
      !end if
      if(isnan(fun))then 
      write(*,*)it,jj
      end if

      
      if(fmax.lt.fun)fmax=fun
      if(fmin.gt.fun)fmin=fun
      fun=fun*wgt/calls
      if(fun.ne.0.)nzero=nzero+1
      fun2=fun*fun
      weight=wgt/calls
c     if(igraph.ne.0)call xbook(0,fun,weight,igraph)
      ti=ti+fun
      sfun2=sfun2+fun2
      do 16 j=1,ndim
      iaj=ia(j)
      d(iaj,j)=d(iaj,j)+fun2
16    continue
11    continue
c.............................................
c  computing the integral and error values
c.............................................
      if (sfun2.ne.0.)go to 1330
      print 1331
      stop
1330  continue
      ti2=ti*ti
      tsi=sqrt((sfun2*calls-ti2)/(calls-1.))
      
      wgt=ti2/tsi**2
      si=si+ti*wgt
      si2=si2+ti2
      swgt=swgt+wgt
      schi=schi+ti2*wgt
      scalls=scalls+calls
      avgi=si/swgt
      yout=avgi
      
      sd=swgt*it/si2
      deint=avgi
      chi2a=0.
      if(it.gt.1)chi2a=sd*(schi/swgt-avgi*avgi)/(it-1)
      sd=1./sqrt(sd)
      err=sd*100./avgi
      yout_err=sd
      yout_sd=chi2a
      write(*,*)it,yout,yout_err,yout_sd
      nzerot=nzerot+nzero
      if(fmaxt.lt.fmax)fmaxt=fmax
      if(fmint.gt.fmin)fmint=fmin
      it0=it
c.............................................
c  printing and histogramming
c.............................................
      if(nprn.ge.2)print 201,it,avgi,sd,err
      if(nprn.ge.3)print 211,ti,tsi,nzero,fmin,fmax,chi2a
      if(nprn.lt.5) go to 21
      do 20 j=1,ndim
      print 202,j
      xin(1)=xi(1,j)
      do 2020 l=2,nd
2020  xin(l)=xi(l,j)-xi(l-1,j)
20    print 204,(xi(i,j),xin(i),d(i,j),i=1,nd)
21    continue
      if(abs(sd/avgi).gt.abs(acc).and.it.lt.itmx)go to 98
      npit=it*ncall
      eff=nzerot*100./npit
      if(nprn.ge.2)print 776,npit,nzerot,eff,fmint,fmaxt
      print 777,avgi,sd,chi2a
c     if(igraph.ne.0)call bookit(2,fun,weight,igraph)
      return
98    continue
c     if(igraph.ne.0)call bookit(0,fun,weight,igraph)
c===========================================================
c c))  redefining the grid
c===========================================================
c.............................................
c  smoothing the f**2 valued stored for each interval
c.............................................
      do 23 j=1,ndim
      xo=d(1,j)
      xn=d(2,j)
      d(1,j)=(xo+xn)/2.
      x(j)=d(1,j)
      do 22 i=2,ndm
      d(i,j)=xo+xn
      xo=xn
      xn=d(i+1,j)
      d(i,j)=(d(i,j)+xn)/3.
      x(j)=x(j)+d(i,j)
22    continue
      d(nd,j)=(xn+xo)/2.
      x(j)=x(j)+d(nd,j)
23    continue
c.............................................
c  computing the 'importance function' of each interval
c.............................................
      do 28 j=1,ndim
      rc=0.
      do 24 i=1,nd
      r(i)=0.
      if(d(i,j).le.0.) go to 224
	IF(D(I,J).LE.1.E-30) THEN
	R(I)=( (X(J)-D(I,J))/X(J)/(LOG(X(J))-LOG(D(I,J))))**ALPH
	ELSE
      xo=x(j)/d(i,j)
      r(i)=((xo-1.)/xo/log(xo))**alph
	END IF
224   rc=rc+r(i)
24    continue
c.............................................
c  redefining the size of each interval
c.............................................
      rc=rc/xnd
      k=0
      xn=0.
      dr=0.
      i=0
25    k=k+1
      dr=dr+r(k)
      xo=xn
      xn=xi(k,j)
26    if(rc.gt.dr) go to 25
      i=i+1
      dr=dr-rc
      xin(i)=xn-(xn-xo)*dr/r(k)
      if(i.lt.ndm) go to 26
      do 27 i=1,ndm
      xi(i,j)=xin(i)
27    continue
      xi(nd,j)=1.
28    continue
c
      go to 9
c===========================================================
c d))  formats for the printouts
c===========================================================
407   format('1  %%%%  routine "vegas" v.1-85 from p.g.lepage',
     . '    (mod.by xd zhang)'/)
290   format('0  %%%%  integration parameters :'/
     . '   # dimensions                =',i8/,
     . '   # calls to f per iteration  =',i8/,
     . '   # iterations maximum        =',i8/,
     . '   # bins in each dimension    =',i8)
201   format(//' iter. no',i3,' acc.results==> int =',g14.5,'+/-',g12.4,
     . '   % error=',g10.2)
211   format(/20x,'iter.results=',g14.5,'+/-',g12.4,
     . '   (f=/=0)=',i6/20x,'   fmin=',g12.4,'   fmax=',g12.4,
     . '   chi**2=',g10.2)
202   format(14h0data for axis,i2 /
     . 7x,'x',9x,'delt x',6x,'sig(f2)',
     . 13x,'x',9x,'delt x',6x,'sig(f2)'/)
204   format(1x,3g12.4,5x,3g12.4)
776   format(//'  %%%% final information : in total '/
     . '  #function calls =',i6,'  #(f=/=0) =',i6,'  % =',g12.4/
     . 5x,'  fmin=',g12.4,'   fmax=',g12.4)
777   format(' '//' ',30('+'),' final result ',30('+')//
     . ' integral value =',g14.5,'+/-',g12.4,6x,
     . ' ( chi**2=',g12.4,')'//' ',74('+'))
1331  format(//
     . '  ## warning: in "vegas" the value of the integral'/
     . '  is exactly zero, so no evolution of the density '/
     . '  distribution can be expected ---> program stop ')

      end

