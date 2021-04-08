% Hua-sheng XIE, huashengxie@gmail.com, CCF-ENN, 2019-06-20 10:13
% Solve Grad-Shafranov equation for FRC equilibrium, with fixed boundary
% condition. Test for Tri alpha's C1.
% 19-08-07 07:07 add finite difference order, =2 (3-point), =4 (5-point)
% 17:50 test ok. To check the accuary, i.e., only 10^-2?
close all;clear all;clc;
irun=0; 
if(irun==0)
  mu0=4*pi*1e-7;
  
  init=0;
  ifd=2; % finite difference order, =2 (3-point), =4 (5-point)
  ibc=2; % b.c. order
  rw =0.17;zw=1*0.75;     %NUCTE-III
  nR=2^7+1; nZ=2^8+1; % grid numbers
  eta=0.5;

  Rlow=0.0*rw; Rup=rw; Zlow=-1.0*zw; Zup=1.0*zw;
  rr=linspace(Rlow,Rup,nR); zz=linspace(Zlow,Zup,nZ);
  dR=rr(2)-rr(1); dZ=zz(2)-zz(1);
  [R,Z]=ndgrid(rr,zz); dR2=dR*dR; dZ2=dZ*dZ;
  Z0=(size(R,2)+1)/2;
  
  load('psib.mat');
  zzb=zz;
  rrb=interp1(zb,rb,zzb,'spline');
  psib=interp1(zb,psib,zzb,'spline');
  
  jrb=round(rrb/dR)+1;
  
  Izeta=-25.0e5;
  psiw=max(abs(psib));
  
  if(init==0)
    C=0.0831875; C1=15623.4955718; C2=4360.3755823;
    Ln=0.7;h=0.7;
%     psi=(0.5-exp(-((R/rw-0.2).^2/0.3^2+((Z+0.0*zw)/zw).^2/0.2^2)))*psiw; % initial psi, to update
  else
    load('psiq1.mat'); % 18-11-13 09:16
%     C=Cq; C1=C1q;C2=C2q;
    C=0.0831875; C1=15623.4955718; C2=4360.3755823;
    Ln=Lnq; Lntemp=Lntemps;
    h=hq;htemp=htemps;
    psi=interp2(Zq,Rq,psiq,Z,R,'nearest');
  end
  
  Ctmp=C; C1tmp=C1; C2tmp=C2;
  
  for jz=1:nZ % fixed boundary condition 
      for jr=jrb(jz):nR
          psi(jr,jz)=psib(jz);
      end
  end
  
  psi0=psi; psi1=psi; psi2=psi;
  dpsim=max(max(abs(psi-psi0)))+0.1;

  it=1; % iterative times
  ita=0;
  runtime0=0;
  
  savepath='./saveit/';
  savefig=0;
  if(savefig==1)
    if(~exist(savepath,'dir')) % in case savepath not exist
      mkdir(savepath);
    end
  end

else
  ita=it;
  runtime0=runtime;
end
runtime=cputime;
nt=70000; % max iterative steps
tol=1e-9; % tolerance for stop
%%
hf=figure('unit','normalized','Position',[0.6 0.35 0.4 0.4],...
   'DefaultAxesFontSize',15);
set(gcf,'PaperPositionMode','auto');

icase =3;
while(it<=nt+ita && dpsim>tol) % iteration, G-S eq
    if icase==1
        dP=-2*C1*C*tanh(C1*(psi-C2)).*sech(C1*(psi-C2)).^2/mu0;
    elseif icase ==2
        dP=-C1*C.*exp(-C1*psi)/mu0;
    elseif icase ==3
        dP=C.*exp(-C1*psi).*(2*C2.*(C2.*psi+1)-C1.*(C2.*psi+1).^2)/mu0;
    end
      
    Jzeta=R.*dP;
    
    % Jztotal=sum(sum(Jzeta))*dr*dz; % not accurate, old version
    Jztotal=0.0;
    for jz=1:nZ % 19-06-23 22:00 update
        for jr=1:(jrb(jz)-1)
            Jztotal=Jztotal+Jzeta(jr,jz)*dR*dZ;
        end
    end
    
    if icase==1
        P=C*sech(C1*(psi-C2))/mu0;  
    elseif icase ==2
        P=C*exp(-C1*psi)/mu0;
    elseif icase ==3
        P=C*exp(-C1*psi).*(C2.*psi+1).^2/mu0;
    end
    
  if(it==1)
    w=0.9; % parameter for iterative, =0.9 for default    
    w0=w;
  else
    wc=w;
    
    P_m=P(:,Z0); psi_m=psi(:,Z0);jr_m=Jzeta(:,Z0);
    lp=solve_lp(P_m,psi_m,rr');
    Lntemp=-(eta+1)*lp;
    fprintf('Ln =%5.4f \n',Lntemp);
    
    htemp=solve_h(jr_m,rr,psi_m);
    fprintf('h =%5.4f \n',htemp);
    
    rs=solve_rs(psi_m,rr);
    fprintf('rs =%5.4f \n',rs);
    
%     C=wc*C+(1-wc)*Ctmp*Izeta/Jztotal;
%     C1=wc*C1+(1-wc)*C1tmp*Lntemp/Ln;
%     C2=wc*C2+(1-wc)*C2tmp*htemp/h;
  end
  
  rhs=psi0;
  for j=2:nZ-1
      for i=2:(jrb(j)-1)
          if(ifd==1) % previous version
              rhs(i,j)=0.5*dZ2*(rr(i)/dR2*((psi0(i+1,j)-...
                  psi0(i,j))/(rr(i)+0.5*dR)-...
                  (psi0(i,j)-psi0(i-1,j))/(rr(i)-0.5*dR))+...
                  (psi0(i,j+1)+psi0(i,j-1))/dZ2+...
                  mu0*rr(i)^2*dP(i,j) ); % 18-11-04 16:25, psi should be psi0
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

  for j=2:nZ-1
    for i=2:(jrb(j)-1)
        psi(i,j)=w*psi0(i,j)+rhs(i,j)*(1-w);
    end
  end
  
  % boundary condition
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
  
  it=it+1;
  fprintf('it = %5.4d \n',it);
  dpsim=max(max(abs(psi-psi0)));
  psi0=psi;
  Ctmp=C; C1tmp=C1; C2tmp=C2;
  
  if(mod(it-2,100)==0) % plot
%%
    ido=find(psi==min(min(psi)));

    subplot(2,3,1:3); hold off;
    %plot(Z,R,'k-',Z',R','k-'); hold on; % 18-11-13 10:43
    indtmp= psi>0; psio=psi; psio(indtmp)=NaN;
    indtmp=find(psi<0); psii=psi; psii(indtmp)=NaN;
    contour(Z,R,psio,10,'linewidth',2); hold on; % insidet %psi>0
    contour(Z,R,psii,15,'linewidth',2); hold on;% outside   %psi<0
    contour(Z,R,psi,[0,0],'r--','linewidth',3); hold on; % highligth separatix 
    plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3); % O-point
    text(0.01*zw+Z(ido(1)),R(ido(1)),'O point','fontsize',13);colorbar;
    
    plot(zb,rb,'-',zzb,rrb,':','linewidth',3);
    
    title(['(a) \Psi(R,Z), FRC G-S equilibrium',...
        ', it=',num2str(it),'/',num2str(nt),...
        ', max(\Psi^{n}-\Psi^{n-1})=',num2str(dpsim)]);
    
    ylabel(['R, r_w=',num2str(rw)]);
    
	strp='p(\psi)=-C*exp(-C1*psi)*(C2*psi+1)^2';
    xlabel(['Z, z_w=',num2str(zw),', I_z=',num2str(Izeta/1e6),...
        'MA, J_{zt}=',num2str(Jztotal/1e6,3),...
        'MA, C=',num2str(C),', C1=',num2str(C1),', C2=',num2str(C2),...
        ' with ',strp]);

%     axis equal;
    box on;
    leg=legend('\psi_{in}','\psi_{out}','separatix');
    legend('boxoff');
    set(leg,'fontsize',8);
    hold off;
    
    subplot(234);
    psi_x=psi(:,Z0);
    ind_min=find(psi_x==min(psi_x));
    plot(rr/rr(ind_min(1)),-Jzeta(:,Z0)/max(abs(Jzeta(:,Z0))),'m','linewidth',3);box on;grid on;
%     plot(rr/rr(ind_min(1)),-dP(:,Z0),'m','linewidth',3);box on;grid on;
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
    xlabel('r'); ylabel('J_{\theta}');

    subplot(235);
    contourf(R,Z,dP); title(['(c) dP, tol=',num2str(tol)]);
    xlabel('R'); ylabel('Z'); box on;

    subplot(236);

    plot(rr/rr(ind_min(1)),psi(:,Z0),'m','linewidth',3);box on;grid on;
    title(['(d) \Psi, w=',num2str(w)]);
    xlabel('r'); ylabel('\psi');
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
    
    pause(0.0001);
    if(savefig==1)
      print(gcf,'-dpng',[savepath,'gs_frc_it=',num2str(it,'%06d'),'.png']);
    end
    %%
    
  end
    
end
runtime=runtime0+cputime-runtime;
fprintf('**************************************** \n');
fprintf('C1 = %7.6f \n',C);
fprintf('C1 = %5.3f \n',C1);
fprintf('C2 = %5.3f \n',C2);

%%
Rq=R;Zq=Z;psiq=psi;C1q=C1;Cq=C;C2q=C2;Lnq=Ln;Lntemps=Lntemp;hq=h;htemps=htemp;
if(Izeta==0)
    filename='psiq_vac';
else
    filename='psiq';
end
save([filename,'.mat'],'Rq','Zq','psiq','Jzeta','P','jrb','zzb','rrb','psib','Cq','C1q',...
    'C2q','nR','nZ','Izeta','hq','htemps');
save([filename,'_ifd=',num2str(ifd),'_ibc=',num2str(ibc),'_nR=',num2str(nR),'_nZ=',num2str(nZ)],...
    'Rq','Zq','psiq','Jzeta','P','jrb','zzb','rrb','psib','Cq','C1q',...
    'C2q','nR','nZ','Izeta');
save('psiq1.mat','Rq','Zq','psiq','Cq','C1q','C2q','Lnq','Lntemps','hq','htemps');
subplot(234);zlabel(['\Psi, runtime=',num2str(runtime),'s']);
set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);
str=['gs_frc_rw=',num2str(rw),',zw=',num2str(zw),...
  ',Iz=',num2str(-Izeta),',nR=',num2str(nR),...
  ',nZ=',num2str(nZ),',C=',num2str(C),',w=',num2str(w),',tol=',num2str(tol),...
  ',it=',num2str(it)];




