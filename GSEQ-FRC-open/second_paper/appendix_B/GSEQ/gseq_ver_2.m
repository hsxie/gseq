% Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, 2019-06-20 10:13
% Solve Grad-Shafranov equation for FRC equilibrium, with fixed boundary
% condition. Test for Tri alpha's C2.
% 19-08-07 07:07 add finite difference order, =2 (3-point), =4 (5-point)
% 17:50 test ok. To check the accuary, i.e., only 10^-2?
close all;

irun=0; % =0, initial run; =1, from previous run
if(irun==0)
  mu0=4*pi*1e-7;
  
  init=0;
  ifd=2; % finite difference order, =2 (3-point), =4 (5-point)
  ibc=2; % b.c. order
  rw=0.17; 
  zw=1.5;
  P0 =1.0/(2*mu0*rw^4);
  
  nR=129; nZ=385; % grid numbers
  
  prs =5.5;     %分界线处压强
  prw =1.0d-5;  %壁压强
  psw =-1.0;    %壁上磁通(外部+等离子)
  pm =-10000.0;  %参数
  e0 =5.8824;
  a0 =prs;
  
  h =1.0/(nR-1);
  fact1 =2.0*(1.0+1.0/e0^2);
  fact2 =0.5*h^2;
  
  a1 =prs/psw*log(prw/prs);
  a2 =0.5*prs*(log(prw/prs)/psw)^2;
  a3 =pm;

  Rlow=0.0*rw; Rup=rw; Zlow=-1.0*zw; Zup=1.0*zw;
  rr=linspace(Rlow,Rup,nR); zz=linspace(Zlow,Zup,nZ);
  dR=rr(2)-rr(1); dZ=zz(2)-zz(1);
  [R,Z] =ndgrid(rr,zz); dR2 =dR*dR; dZ2 =dZ*dZ;
  
  load('psib.mat');
  zzb =zz;
  rrb =interp1(zb,rb,zzb,'spline');
  psib =interp1(zb,psib,zzb,'spline');
  
  jrb =round(rrb/dR)+1;
  
  Izeta =-6.0e5;
%   Izeta=-0.0e5;
  psiw =max(abs(psib));
  
  if(init==0)
    C=0.0067; C3=-0.0005;
    psi=0.0;
%     psi=-(0.5-exp(-((R/rw-0.25).^2/0.2^2+(Z/zw).^2/0.2^2)))*psiw; % initial psi, to update
  else
    load('psiq1.mat'); % 18-11-13 09:16
    C=Cq; C3=C3q;
    psi=interp2(Zq,Rq,psiq,Z,R,'nearest');
  end

  Ctmp=C; Ctmp3=C3;

% bound condition
for jz=1:nZ
    psi(nR,jz)=psib(jz);
end
  
  psi0=psi; psi1=psi; psi2=psi;
  dpsim=max(max(abs(psi-psi0)))+0.1;

  it=1; % iterative times
  ita=0;
  runtime0=0;

else
  ita=it;
  runtime0=runtime;
%     runtime0=150000;
end
runtime=cputime;
nt=150000; % max iterative steps
tol=1e-10; % tolerance for stop
%%
hf=figure('unit','normalized','Position',[0.5 0.1 0.5 0.4],...
   'DefaultAxesFontSize',15);
while(it<=nt+ita && dpsim>tol) % iteration, G-S eq
   
    psi0=psi;
    C2=2000;
    mcase =2;
    if mcase ==1
        dP =-C2*C*tanh(C2*(psi-C3)).*sech(C2*(psi-C3))/mu0;
        P=C*sech(C2*(psi-C3))/mu0;
    elseif mcase ==2
        for i=1:nR
            for j=1:nZ
                if (psi(i,j)>=0.0)
                    dP(i,j) =(a1+2.0*a2*psi(i,j)+3.0*a3*psi(i,j)*psi(i,j));  
                    P(i,j) =(a0+a1*psi(i,j)+a2*psi(i,j)^2+a3*psi(i,j)^3);
                else 
                    dP(i,j) =(prs/psw*(prw/prs)^(psi(i,j)/psw)*log(prw/prs));
                     P(i,j) =(prs*(prw/prs)^(psi(i,j)/psw));
                end
            end
        end
    end
    
    Jzeta=R.*dP;
    Jztotal=0.0;
    
    for jz=1:nZ % 19-06-23 22:00 update
        for jr=1:(jrb(jz)-1)
%          for jr=1:nR
            Jztotal=Jztotal+Jzeta(jr,jz)*dR*dZ;
        end
    end
    
  if(it==1)
    w=0.5;
    w0=w;
    om =1.95;
  else
    wc=w;
    C=wc*C+(1-wc)*Ctmp*Izeta/Jztotal;
    om =1.95;
  end
  
  rhs=psi0;
  cal_case =2; 
  if cal_case ==1
      for j=2:nZ-1
          for i=2:(jrb(j)-1)
              if(ifd==1) % previous version
                  rhs(i,j)=0.5*dZ2*( rr(i)/dR2*((psi0(i+1,j)-...
                      psi0(i,j))/(rr(i)+0.5*dR)-...
                      (psi0(i,j)-psi0(i-1,j))/(rr(i)-0.5*dR))+...
                      (psi0(i,j+1)+psi0(i,j-1))/dZ2+...
                      fact2*mu0*rr(i)^2*dP(i,j) ); % 18-11-04 16:25, psi should be psi0
              elseif(ifd==2) % 3-point
                  rhs(i,j)=1/(2/dR2+2/dZ2)*(...
                      (psi0(i+1,j)+psi0(i-1,j))/dR2-...
                      (psi0(i+1,j)-psi0(i-1,j))/(2*dR*rr(i))+...
                      (psi0(i,j+1)+psi0(i,j-1))/dZ2+...
                      mu0*rr(i)^2*dP(i,j));
              elseif(ifd==4) % 5-point
                  if((j>2 && j<nZ-1) && (i>2 && i<jrb(j)-1))
                    rhs(i,j)=1/(30/(12*dR2)+30/(12*dZ2))*(...
                      (-psi0(i+2,j)+16*psi0(i+1,j)+16*psi0(i-1,j)-psi0(i-2,j))/(12*dR2)-...
                      (-psi0(i+2,j)+8*psi0(i+1,j)-8*psi0(i-1,j)+psi0(i-2,j))/(12*dR*rr(i))+...
                      (-psi0(i,j+2)+16*psi0(i,j+1)+16*psi0(i,j-1)-psi0(i,j-2))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  elseif((j>2 && j<nZ-1) && (i==2))
                    rhs(i,j)=1/(20/(12*dR2)-10/(12*dR*rr(i))+30/(12*dZ2))*(...
                      (-psi0(i+3,j)+4*psi0(i+2,j)+6*psi0(i+1,j)+11*psi0(i-1,j))/(12*dR2)-...
                      (psi0(i+3,j)-6*psi0(i+2,j)+18*psi0(i+1,j)-3*psi0(i-1,j))/(12*dR*rr(i))+...
                      (-psi0(i,j+2)+16*psi0(i,j+1)+16*psi0(i,j-1)-psi0(i,j-2))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  elseif((j==2) && (i==2))%
                    rhs(i,j)=1/(20/(12*dR2)-10/(12*dR*rr(i))+20/(12*dZ2))*(...
                      (-psi0(i+3,j)+4*psi0(i+2,j)+6*psi0(i+1,j)+11*psi0(i-1,j))/(12*dR2)-...
                      (psi0(i+3,j)-6*psi0(i+2,j)+18*psi0(i+1,j)-3*psi0(i-1,j))/(12*dR*rr(i))+...
                      (-psi0(i,j+3)+4*psi0(i,j+2)+6*psi0(i,j+1)+11*psi0(i,j-1))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  elseif((j==nZ-1) && (i==2))%
                    rhs(i,j)=1/(20/(12*dR2)-10/(12*dR*rr(i))+20/(12*dZ2))*(...
                      (-psi0(i+3,j)+4*psi0(i+2,j)+6*psi0(i+1,j)+11*psi0(i-1,j))/(12*dR2)-...
                      (psi0(i+3,j)-6*psi0(i+2,j)+18*psi0(i+1,j)-3*psi0(i-1,j))/(12*dR*rr(i))+...
                      (-psi0(i,j-3)+4*psi0(i,j-2)+6*psi0(i,j-1)+11*psi0(i,j+1))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  elseif((j>2 && j<nZ-1) && (i==(jrb(j)-1)))
                    rhs(i,j)=1/(20/(12*dR2)+10/(12*dR*rr(i))+30/(12*dZ2))*(...
                      (-psi0(i-3,j)+4*psi0(i-2,j)+6*psi0(i-1,j)+11*psi0(i+1,j))/(12*dR2)-...
                      (-psi0(i-3,j)+6*psi0(i-2,j)-18*psi0(i-1,j)+3*psi0(i+1,j))/(12*dR*rr(i))+...
                      (-psi0(i,j+2)+16*psi0(i,j+1)+16*psi0(i,j-1)-psi0(i,j-2))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  elseif((j==2) && (i==(jrb(j)-1)))%
                    rhs(i,j)=1/(20/(12*dR2)+10/(12*dR*rr(i))+20/(12*dZ2))*(...
                      (-psi0(i-3,j)+4*psi0(i-2,j)+6*psi0(i-1,j)+11*psi0(i+1,j))/(12*dR2)-...
                      (-psi0(i-3,j)+6*psi0(i-2,j)-18*psi0(i-1,j)+3*psi0(i+1,j))/(12*dR*rr(i))+...
                      (-psi0(i,j+3)+4*psi0(i,j+2)+6*psi0(i,j+1)+11*psi0(i,j-1))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  elseif((j==nZ-1) && (i==(jrb(j)-1)))%
                    rhs(i,j)=1/(20/(12*dR2)+10/(12*dR*rr(i))+20/(12*dZ2))*(...
                      (-psi0(i-3,j)+4*psi0(i-2,j)+6*psi0(i-1,j)+11*psi0(i+1,j))/(12*dR2)-...
                      (-psi0(i-3,j)+6*psi0(i-2,j)-18*psi0(i-1,j)+3*psi0(i+1,j))/(12*dR*rr(i))+...
                      (-psi0(i,j-3)+4*psi0(i,j-2)+6*psi0(i,j-1)+11*psi0(i,j+1))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  elseif((j==2) && (i>2 && i<jrb(j)-1))
                    rhs(i,j)=1/(30/(12*dR2)+20/(12*dZ2))*(...
                      (-psi0(i+2,j)+16*psi0(i+1,j)+16*psi0(i-1,j)-psi0(i-2,j))/(12*dR2)-...
                      (-psi0(i+2,j)+8*psi0(i+1,j)-8*psi0(i-1,j)+psi0(i-2,j))/(12*dR*rr(i))+...
                      (-psi0(i,j+3)+4*psi0(i,j+2)+6*psi0(i,j+1)+11*psi0(i,j-1))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  elseif((j==(nZ-1)) && (i>2 && i<jrb(j)-1))
                    rhs(i,j)=1/(30/(12*dR2)+20/(12*dZ2))*(...
                      (-psi0(i+2,j)+16*psi0(i+1,j)+16*psi0(i-1,j)-psi0(i-2,j))/(12*dR2)-...
                      (-psi0(i+2,j)+8*psi0(i+1,j)-8*psi0(i-1,j)+psi0(i-2,j))/(12*dR*rr(i))+...
                      (-psi0(i,j-3)+4*psi0(i,j-2)+6*psi0(i,j-1)+11*psi0(i,j+1))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));   
                  end
              elseif(ifd==3) % 5-point + 3-point
                  if((j>2 && j<nZ-1) && (i>2 && i<jrb(j)-1))
                    rhs(i,j)=1/(30/(12*dR2)+30/(12*dZ2))*(...
                      (-psi0(i+2,j)+16*psi0(i+1,j)+16*psi0(i-1,j)-psi0(i-2,j))/(12*dR2)-...
                      (-psi0(i+2,j)+8*psi0(i+1,j)-8*psi0(i-1,j)+psi0(i-2,j))/(12*dR*rr(i))+...
                      (-psi0(i,j+2)+16*psi0(i,j+1)+16*psi0(i,j-1)-psi0(i,j-2))/(12*dZ2)+...
                      mu0*rr(i)^2*dP(i,j));
                  else
                      rhs(i,j)=1/(2/dR2+2/dZ2)*(...
                          (psi0(i+1,j)+psi0(i-1,j))/dR2-...
                          (psi0(i+1,j)-psi0(i-1,j))/(2*dR*rr(i))+...
                          (psi0(i,j+1)+psi0(i,j-1))/dZ2+...
                          mu0*rr(i)^2*dP(i,j));
                  end
              end
          end
      end
  elseif cal_case ==2
      for i=2:nR-1
          for j=1:nZ
            if j==1
                rhs(i,j) =1.0/fact1*(psi0(i+1,j)+psi0(i-1,j)-0.5*h/rr(i)*(psi0(i+1,j)-psi0(i-1,j)) ...
                    +1.0/e0^2*(psi0(i,j+1)+psi0(i,j+1))+fact2*rr(i)^2*dP(i,j));
            elseif j==nZ
                rhs(i,j) =1.0/fact1*(psi0(i+1,j)+psi0(i-1,j)-0.5*h/rr(i)*(psi0(i+1,j)-psi0(i-1,j)) ...
                    +1.0/e0^2*(psi0(i,j-1)+psi0(i,j-1))+fact2*rr(i)^2*dP(i,j));
            else
                rhs(i,j) =1.0/fact1*(psi0(i+1,j)+psi0(i-1,j)-0.5*h/rr(i)*(psi0(i+1,j)-psi0(i-1,j)) ...
                    +1.0/e0^2*(psi0(i,j+1)+psi0(i,j-1))+fact2*rr(i)^2*dP(i,j));
            end
          end
      end
  end    
  
  if cal_case ==1
      for j=2:nZ-1
        for i=2:(jrb(j)-1)
            psi(i,j)=w*psi0(i,j)+rhs(i,j)*(1-w);
        end
       end
  elseif cal_case ==2
      for j=1:nZ
        for i=1:nR
            psi(i,j) = psi(i,j)+w*(rhs(i,j)-psi(i,j));
        end
      end
  end
  
  %% boundary condition
  if cal_case ==1
      psi(1,:)=0.0; % zero boundary at axis
      for jz=1:nZ % fixed boundary condition at wall
          psi(jrb(jz),jz)=psib(jz);
      end
  
  % dpsi/dz=0 at z=min,max
      if(ibc==1) % 2-point
          psi(:,1)=psi(:,2);
          psi(:,nZ)=psi(:,nZ-1);
      elseif(ibc==2) % 3-point
          psi(:,1)=(4*psi(:,2)-psi(:,3))/3;
          psi(:,nZ)=(4*psi(:,nZ-1)-psi(:,nZ-2))/3;
      elseif(ibc==3) % 4-point
          psi(:,1)=(18*psi(:,2)-9*psi(:,3)+2*psi(:,4))/11;
          psi(:,nZ)=(18*psi(:,nZ-1)-9*psi(:,nZ-2)+2*psi(:,nZ-3))/11;
      elseif(ibc==4) % 5-point
          psi(:,1)=(48*psi(:,2)-36*psi(:,3)+16*psi(:,4)-3*psi(:,5))/25;
          psi(:,nZ)=(48*psi(:,nZ-1)-36*psi(:,nZ-2)+16*psi(:,nZ-3)-3*psi(:,nZ-4))/25;
      end
  end
  
  it=it+1;
  dpsim=max(max(abs(psi-psi0)));
  psi0=psi;
  Ctmp=C;

  if(mod(it-2,100)==0) % plot
%%
    ido=find(psi==min(min(psi)));

    subplot(2,3,1:3); hold off;
    indtmp= psi>0; psio=psi; psio(indtmp)=NaN;
    indtmp=find(psi<0); psii=psi; psii(indtmp)=NaN;
    contour(Z,R,psio,10,'linewidth',2); hold on;
    contour(Z,R,psii,20,'linewidth',2); hold on;
    contour(Z,R,psi,[0,0],'r--','linewidth',3); hold on; % highligth separatix 
    plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3); % O-point
    text(0.01*zw+Z(ido(1)),R(ido(1)),'O point','fontsize',13);
    colorbar;
    
    title(['(a) \Psi(R,Z), FRC G-S equilibrium',...
        ', it=',num2str(it),'/',num2str(nt),...
        ', max(\Psi^{n}-\Psi^{n-1})=',num2str(dpsim)]);
    
    ylabel(['R, r_w=',num2str(rw)]);
    
	strp='p(\psi)=C*sech(C2*(psi-C3))';
    xlabel(['Z, z_w=',num2str(zw),', I_z=',num2str(Izeta/1e6),...
        'MA, J_{zt}=',num2str(Jztotal/1e6,3),...
        'MA, C=',num2str(C),', C2=',num2str(C2),', C3=',num2str(C3),...
        ' with ',strp]);
    box on;grid on;
    leg=legend('\psi_{in}','\psi_{out}','separatix');
    legend('boxoff');
    set(leg,'fontsize',8);
    hold off;
    
    subplot(2,3,4:5);
    contourf(Z,R,P); title(['(c) P, tol=',num2str(tol)]);
    xlabel('R'); ylabel('Z'); box on;grid on;colorbar;
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);

    subplot(236);
    v=[0,-0.1,0.1,-0.01,0.01];
    contour(Z,R,psi,v); box on;
    title(['(d) \Psi, w=',num2str(w)]);
    xlabel('R'); ylabel('Z');
    
    pause(0.0001);
%     if(savefig==1)
%       print(gcf,'-dpng',[savepath,'gs_frc_it=',num2str(it,'%06d'),'.png']);
%     end
  end
    
end
runtime=runtime0+cputime-runtime;
psiq=psi;
filename='psi';
save([filename,'.mat'],'psiq');
