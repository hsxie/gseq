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
    irs=init;
    icase=1;
    if icase==1
        load('MRR1d.mat');
    elseif icase==2
        load('gs1d.mat');
        alpha=alphas;
        sigma=sigmas;
        n=2;
    end
    
    zw=0.75;
    nR=2^9+1; nZ=2^8+1; % grid numbers
    Be=1.0;   
    
    Pm0=Be^2/(2*mu0);
    
    Rlow=0.0*Rw; Rup=Rw; Zlow=-1.0*zw; Zup=1.0*zw;
    rr=linspace(Rlow,Rup,nR); zz=linspace(Zlow,Zup,nZ);
    dR=rr(2)-rr(1); dZ=zz(2)-zz(1);
    [R,Z]=ndgrid(rr,zz); dR2=dR*dR; dZ2=dZ*dZ;
    Z0=(size(R,2)+1)/2;
    
    load('psib_paper.mat');
    zzb=zz;
    rrb=interp1(zb,rb,zzb,'spline');
    %   psib=interp1(zb,psib.*(Be*Rs^2),zzb,'spline');
    psib=interp1(zb,psib,zzb,'spline');
    
    jrb=round(rrb/dR)+1;
    
    psiw=max(abs(psib));
    drs=1;
    
    if(init==0)
        psi=zeros(nR,nZ);
        for i=1,nR
            psi(i,:)=psi1d(i)*Be*Rs^2;
        end
        C=1;
    else   
        load('psiq.mat');
        psi=interp2(Zq,Rq,psiq,Z,R,'nearest');
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
nt=50000*2; % max iterative steps
tol=1e-13; % tolerance for stop
%%
hf=figure('unit','normalized','Position',[0.6 0.35 0.4 0.4],...
    'DefaultAxesFontSize',15);
while(it<=nt+ita && dpsim>tol)
    if (1==1)
        dP=C*Pm0*betas.*exp(-alpha.*psi./(Be*Rs^2)).*(sigma.*psi./(Be*Rs^2)+1).^(n-1).*...
            (-alpha/(Be*Rs^2).*(sigma.*psi./(Be*Rs^2)+1)+n*sigma/(Be*Rs^2));
    elseif (1==2)
        Rm=C*R0;
        dP=Pm0*betas.*exp(-4*K/(Be*Rm^2).*psi).*(2*4*sigma/(Be*Rm^2).*(4*sigma/(Be*Rm^2).*psi+1)...
            -4*K/(Be*Rm^2).*(4*sigma/(Be*Rm^2).*psi+1).^2);
    end
    
    Jzeta=-R.*dP;
    Jztotal=0.0;
    Jztotal1=0.0;Jztotal2=0;
%     for jz=1:nZ % 19-06-23 22:00 update
    for jz=floor(0.5*nZ) % 19-06-23 22:00 update
        for jr=1:(jrb(jz)-1)
%             Jztotal1=Jztotal1+2*pi*Jzeta(jr,jz)*dR*dZ;
            Jztotal=Jztotal+2*pi*Jzeta(jr,jz)*R(jr,jz)*dR;
        end
    end
    
    for jz=1:nZ % 19-06-23 22:00 update
        for jr=1:(jrb(jz)-1)
            Jztotal1=Jztotal1+Jzeta(jr,jz)*dR*dZ;
            Jztotal2=Jztotal2+2*pi*Jzeta(jr,jz)*R(jr,jz)*dR*dZ;
        end
    end
    
    
    if (1==1)
        P=C*Pm0*betas.*exp(-alpha.*psi./(Be*Rs^2)).*(sigma.*psi./(Be*Rs^2)+1).^n;
    elseif (1==2)
        P=Pm0*betas.*exp(-4*K/(Be*Rm^2).*psi).*(4*sigma/(Be*Rm^2).*psi+1).^2;
    end
    Jzetaz0=Jzeta(:,floor(0.5*nZ));
    Prz0=P(:,floor(0.5*nZ));
    
    if(it==1)
        w=0.1; % parameter for iterative, =0.9 for default
        w0=w;
    elseif(it>=100)
        wc=0.1;
        
%         Izeta=pi*sum(2*Jzeta.*rr*dR);
        Izeta=Jztotal;
        C=C*Izeta0/Izeta;
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
    
    psimid=psi(:,Z0);
    
    if irs==2
        rs=solve_rs(psimid,rr);
        drs=abs(rs-Rs)/Rs;
    elseif (irs==0 || irs==1)
        itmp=floor(0.15*nR);
        rs=rr(find(abs(psimid(itmp:end))==min(abs(psimid(itmp:end))))+itmp-1);
        rs=rs(1);
    end
    
    it=it+1;
    dpsim=max(max(abs(psi-psi0)));
    psi0=psi;
    
    if(mod(it-2,100)==0) % plot
        fprintf('it = %5.4d \n',it);
        fprintf('rs = %5.4d \n',rs);
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
        text(0.01*zw+Z(ido(1)),R(ido(1)),'O point','fontsize',13);colorbar;colormap(hot);
        
        plot(zb,rb,'-',zzb,rrb,':','linewidth',3);
        
        title(['(a) \Psi(R,Z), FRC G-S equilibrium',...
            ', it=',num2str(it),'/',num2str(nt),...
            ', max(\Psi^{n}-\Psi^{n-1})=',num2str(dpsim)]);
        
        ylabel(['R, r_w=',num2str(Rw)]);
        
        strp='p(\psi)=\beta_s exp(\alpha \psi)*(\sigma \psi+1)^2';
        xlabel(['Z, z_w=',num2str(zw),'; J_{zt}=',num2str(Jztotal/1e6,3),'MA; \beta_s =',num2str(betas),...
            '; \alpha =',num2str(alpha),'; \sigma =',num2str(sigma),...
            '; with ',strp]);
        
        leg=legend('\psi_{in}','\psi_{out}','separatix');
        legend('boxoff');set(leg,'fontsize',8);
        set(gca,'xtick',[-0.75:0.25:0.75]);
        hold off;
        
        subplot(234);
        psi_x=psi(:,Z0);
        ind_min=find(psi_x==min(psi_x));
        %     plot(rr/rr(ind_min(1)),Jzeta(:,Z0)/max(abs(Jzeta(:,Z0))),'m','linewidth',3);box on;grid on;
        plot(rr,Jzeta(:,Z0)/max(abs(Jzeta(:,Z0))),'m','linewidth',3);box on;grid on;
        %     plot(rr/rr(ind_min(1)),-dP(:,Z0),'m','linewidth',3);box on;grid on;
        %     xlim([0 2]);
        set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
        xlabel('r/R_0'); ylabel('J_{\theta}');
        
        subplot(235);
        contourf(R,Z,dP); title(['(c) dP, tol=',num2str(tol)]);
        xlabel('R'); ylabel('Z'); box on;
        
        subplot(236);
        
        %     plot(rr/rr(ind_min(1)),psi(:,Z0),'m','linewidth',3);box on;grid on;
        plot(rr,psi(:,Z0),'m','linewidth',3);box on;grid on;
        title(['(d) \Psi, w=',num2str(w)]);
        xlabel('r/R_0'); ylabel('\psi');
        %     xlim([0,1.6]);
        set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
        
        pause(0.0001);
%         if(savefig==1)
%             print(gcf,'-dpng',[savepath,'gs_frc_it=',num2str(it,'%06d'),'.png']);
%         end
        %%
        
    end
    
end
runtime=runtime0+cputime-runtime;

%%
Rq=R;Zq=Z;psiq=psi;
filename='psiq';
save([filename,'.mat'],'Rq','Zq','psiq','Jzeta','P','jrb','zzb','rrb','psib','Z0','Be','Rs',...
    'nR','nZ','C');
save('psiq1.mat','Rq','Zq','psiq','drs','C');
