% Hua-sheng Xie, huashengxie@gmail.com, 2019-05-29 14:47
% Read psi(r,z) and plot
close all; clear; clc;

% load('psiq.mat');
load('psiq.mat')
% load('psiq_I-1.0_ln-040.mat');
rw=max(max(Rq)); zw=max(max(Zq)); zwa=min(min(Zq));
[nRq,nZq]=size(psiq);

nR=nRq; nZ=nZq;
% nR=2^6+1; nZ=2^8+1;
Rlow=0; Rup=rw; Zlow=zwa; Zup=zw;
rr=linspace(Rlow,Rup,nR); zz=linspace(Zlow,Zup,nZ);
dR=rr(2)-rr(1); dZ=zz(2)-zz(1);
[R,Z]=ndgrid(rr,zz); dR2=dR*dR; dZ2=dZ*dZ;
psi=interp2(Zq,Rq,psiq,Z,R,'spline'); % 'cubic'
% psi=interp2(Zq,Rq,psiq,Z,R); % not ok for J&dp, due to 'linear' interp
Jztotal=sum(sum(Jzeta))*dR*dZ;

data_psi=psi<0; Jzx=Jzeta; Jzx(data_psi)=0;
Jztotal_inside=sum(sum(Jzx))*dR*dZ


Br=0.*psi; Bz=0.*psi; dp=0.*psi; J=0.*psi;
for jZ=2:(nZ-1)
    for jR=2:(nR-1)
        Br(jR,jZ)=-1/R(jR,jZ)*(psi(jR,jZ+1)-psi(jR,jZ-1))/(2*dZ);
        Bz(jR,jZ)= 1/R(jR,jZ)*(psi(jR+1,jZ)-psi(jR-1,jZ))/(2*dR);
        J(jR,jZ)=-1/R(jR,jZ)*((psi(jR+1,jZ)-2*psi(jR,jZ)+psi(jR-1,jZ))/(dR^2)...
            -1/R(jR,jZ)*(psi(jR+1,jZ)-psi(jR-1,jZ))/(2*dR)...
            +(psi(jR,jZ+1)-2*psi(jR,jZ)+psi(jR,jZ-1))/(dZ^2)); % to check
%         J(jR,jZ)=((psi(jR,jZ+1)-2*psi(jR,jZ)+psi(jR,jZ-1))); % to check
        dp(jR,jZ)=J(jR,jZ)/R(jR,jZ);
    end
end

for jZ=1:nZ
    psi(1,jZ)=psi(2,jZ); psi(nR,jZ)=psi(nR-1,jZ);
    Br(1,jZ)=Br(2,jZ); Br(nR,jZ)=Br(nR-1,jZ);
    Bz(1,jZ)=Bz(2,jZ); Bz(nR,jZ)=Bz(nR-1,jZ);
    J(1,jZ)=J(2,jZ); J(nR,jZ)=J(nR-1,jZ);
    dp(1,jZ)=dp(2,jZ); dp(nR,jZ)=dp(nR-1,jZ);
    
    for jR=(jrb(jZ)+1):nR
%         psi(jR,jZ)=0.0;
%         Br(jR,jZ)=0.0;
%         Bz(jR,jZ)=0.0;
%         J(jR,jZ)=0.0;
%         dp(jR,jZ)=0.0;
        psi(jR,jZ)=NaN;
        Br(jR-1,jZ)=NaN;
        Bz(jR-1,jZ)=NaN;
        J(jR-1,jZ)=NaN;
        dp(jR-1,jZ)=NaN;
    end
end
% dp(dp>1e3)=NaN;

[jRo,jZo]=find(psi==min(min(psi)));
[jRm,jZm]=find(psi==max(max(psi)));
ro=R(jRo(1),jZo(1)); zo=Z(jRo(1),jZo(1));
zm=Z(jRm(1),jZm(1));
% rm=R(jRm(1),jZm(1));
% rm=0.14;
rm=0.15;

%%
close all;
hf=figure('unit','normalized','Position',[0.6 0.1 0.4 0.7],...
   'DefaultAxesFontSize',12);
hs1=subplot(611);
% set(hs1,'position',[0.1,0.7,0.8,0.1]);
pcolor(Z,R,psi); shading interp; hold on; colorbar('northOutside');
set(gca,'Fontsize',17);
% vv=[0,0.01,-0.001,-0.01,-0.005,0.15,0.4];
psimin=min(min(psi));
vv=[0,0.01,0.15,0.4,0.9*psimin,0.3*psimin,0.03*psimin];
[C1,h1]=contour(Z,R,psi,vv,'linewidth',1,'linestyle','-',...
    'color','k','visible','off'); hold on;
[C2,h2]=contour(Z,R,psi,[0,0],'linewidth',1,'linestyle',...
    '--','color','r','visible','off'); hold on;
ylabel('R (m)'); axis equal; axis tight; box on; % title('(a) \psi');
plot(zo,ro,'rx');
ylim([0,rm]);
set(gca,'ytick',[0:0.2:rm]);

% grid on;
subplot(612);
pcolor(Z,R,Br); shading interp; hold on; colorbar('northOutside');
set(gca,'Fontsize',17);
ylabel('R (m)'); axis equal; axis tight; box on; % title('(b) B_r');
ylim([0,rm]);
set(gca,'ytick',[0:0.2:rm]);
subplot(613);
pcolor(Z,R,Bz); shading interp; hold on; colorbar('northOutside');
set(gca,'Fontsize',17);
ylabel('R (m)'); axis equal; axis tight; box on; % title('(c) B_z');
set(gca,'ytick',[0:0.2:rm]);

subplot(614);
% contour(Z,R,J);
pcolor(Z,R,J); shading interp; hold on; colorbar('northOutside');
set(gca,'Fontsize',17);
ylabel('R (m)'); axis equal; axis tight; box on; % title('(d) J_\zeta');
set(gca,'ytick',[0:0.2:rm]);
subplot(615);
% dp(abs(dp)>2e1)=NaN;
pcolor(Z,R,dp); shading interp; hold on; colorbar('northOutside');
set(gca,'Fontsize',17);
ylabel('R (m)'); axis equal; axis tight; box on;%title('(e) dp/d\psi');
set(gca,'ytick',[0:0.2:rm]);
% axis off;

JBz=J.*Bz;
JBr=J.*Br;
%
% subplot(511);
% hold on;plot(zo,ro,'x',zm,rm,'+'); hold on;

% calculate p(r,z), 19-05-30 15:36
imethod=2;
p=0.*R;
if(imethod==1) % calculate fro dpsi

    for jR=2:nR
        p(jR,jZo)=p(jR-1,jZo)+0.5*(dp(jR-1,jZo)+dp(jR,jZo))*(psi(jR,jZo)-psi(jR-1,jZo));
    end
    for jR=1:nR
        for jZ=(jZo-1):-1:1
            p(jR,jZ)=p(jR,jZ+1)+0.5*(dp(jR,jZ)+dp(jR,jZ+1))*(psi(jR,jZ)-psi(jR,jZ+1));
        end
        for jZ=(jZo+1):1:nZ
            p(jR,jZ)=p(jR,jZ-1)+0.5*(dp(jR,jZ)+dp(jR,jZ-1))*(psi(jR,jZ)-psi(jR,jZ-1));
        end
    end
else % calculate from dr and dz
    for jR=2:nR
        p(jR,jZo)=p(jR-1,jZo)+0.5*(JBz(jR-1,jZo)+JBz(jR,jZo))*dR;
    end
    for jR=1:nR
        for jZ=(jZo-1):-1:1
            p(jR,jZ)=p(jR,jZ+1)+0.5*(JBr(jR,jZ)+JBr(jR,jZ+1))*dZ;
        end
        for jZ=(jZo+1):1:nZ
            p(jR,jZ)=p(jR,jZ-1)-0.5*(JBr(jR,jZ)+JBr(jR,jZ-1))*dZ;
        end
    end
end
p=p-min(min(p))+0.05*(max(max(p)-min(min(p))));
subplot(6,1,6);

pcolor(Z,R,p); shading interp; hold on; colorbar('northOutside');
set(gca,'Fontsize',17);
xlabel('Z (m)');ylabel('R (m)');axis equal; axis tight; box on;%title('(f) p(\psi)'); 
set(gca,'ytick',[0:0.2:rm]);
% contour(Z,R,p);axis equal;

% axis off;


%%

C1(C1>5)=NaN;
C2(C2>5)=NaN;
% strj={'\psi','B_r','B_z','J_\zeta','dp/d\psi','p(\psi)'};
strj={'\psi','B_r','B_z','J_\zeta','dp/d\psi','p(\psi)'};
for j=1:6
    subplot(6,1,j);
    plot(C1(1,:),C1(2,:),'k','linewidth',1);
    plot(C2(1,:),C2(2,:),'r--','linewidth',1);
    stj=['(',char(96+j),') '];
%     stj=strcat(stj,' ',strj(j));
    stj=strcat(strj(j));
%     text(zw,-0.15*rw,stj,'fontsize',14);
    text(0.8*zw,0.8*rw,stj,'fontsize',18,'color','m');
end

% set(gcf,'PaperPositionMode','auto');
% print(gcf,'-dpng','equilibrium.png');
%%
dpdr=0.*p;
dpdz=0.*p;
for jR=2:(nR-1)
    for jZ=2:(nZ-1)
        dpdr(jR,jZ)=(p(jR+1,jZ)-p(jR-1,jZ))/(2*dR);
        dpdz(jR,jZ)=(p(jR,jZ+1)-p(jR,jZ-1))/(2*dZ);
    end
end
JBr_o_dpdz=-JBr./dpdz;
JBz_o_dpdr=JBz./dpdr;

[jr,jz]=ndgrid(1:nR,1:nZ);
eqdat=[reshape(jr,1,[]);reshape(jz,1,[]);reshape(R,1,[]);...
    reshape(Z,1,[]);reshape(psi,1,[]);reshape(p,1,[]);reshape(J,1,[])];
% fid = fopen('frceq.dat', 'w');
% fprintf(fid,'%d\n%d\n',nR,nZ);
% fprintf(fid,'%5d %5d %16.7e %16.7e %16.7e %16.7e %16.7e \n',eqdat);
% fclose(fid);
