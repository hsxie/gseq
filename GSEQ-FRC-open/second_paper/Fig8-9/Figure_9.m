close all;clear;clc;
fontsize=25;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.45,0.4,0.5,0.4]);

%%%%%%%%%%%%%%%%%%%%%%%deltas=0.6%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('psiq1_fixed_betas=0_4_deltas=0_15.mat');

rw=max(max(Rq)); zw=max(max(Zq)); zwa=min(min(Zq));
[nRq,nZq]=size(psiq);

nR=nRq; nZ=nZq;
Rlow=0; Rup=rw; Zlow=zwa; Zup=zw;
rr=linspace(Rlow,Rup,nR); zz=linspace(Zlow,Zup,nZ);
dR=rr(2)-rr(1); dZ=zz(2)-zz(1);
[R,Z]=ndgrid(rr,zz); dR2=dR*dR; dZ2=dZ*dZ;
psi_betas4=interp2(Zq,Rq,psiq,Z,R,'spline'); % 'cubic'

for jZ=1:nZ
    psi_betas4(1,jZ)=psi_betas4(2,jZ); psi_betas4(nR,jZ)=psi_betas4(nR-1,jZ);
end
Pr_betas4=P(:,Z0)./max(P(:,Z0));

%%%%%%%%%%%%%%%%%%%%%%%deltas=0.7%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('psiq1_fixed_betas=0_5_deltas=0_15.mat');
psi_betas5=interp2(Zq,Rq,psiq,Z,R,'spline'); % 'cubic'
for jZ=1:nZ
    psi_betas5(1,jZ)=psi_betas5(2,jZ); psi_betas5(nR,jZ)=psi_betas5(nR-1,jZ);
end
Pr_betas5=P(:,Z0)./max(P(:,Z0));
%%%%%%%%%%%%%%%%%%%%%%%deltas=0.8%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('psiq1_fixed_betas=0_6_deltas=0_15.mat');
psi_betas6=interp2(Zq,Rq,psiq,Z,R,'spline'); % 'cubic'
for jZ=1:nZ
    psi_betas6(1,jZ)=psi_betas6(2,jZ); psi_betas6(nR,jZ)=psi_betas6(nR-1,jZ);
end
Pr_betas6=P(:,Z0)./max(P(:,Z0));
%%%%%%%%%%%%%%%%%%%%%%%deltas=0.9%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('psiq1_fixed_betas=0_7_deltas=0_15.mat');
psi_betas7=interp2(Zq,Rq,psiq,Z,R,'spline'); % 'cubic'
for jZ=1:nZ
    psi_betas7(1,jZ)=psi_betas7(2,jZ); psi_betas7(nR,jZ)=psi_betas7(nR-1,jZ);
end
Pr_betas7=P(:,Z0)./max(P(:,Z0));

axes('position',[0.065,0.58,0.91,0.4]);%fig-1[left bottom width height]
ido=find(psi_betas4==min(min(psi_betas4)));
indtmp= psi_betas4>0; psio=psi_betas4; psio(indtmp)=NaN;
indtmp=find(psi_betas4<0); psii=psi_betas4; psii(indtmp)=NaN;
contour(Z,R,psio,12,'linewidth',2); hold on;
contour(Z,R,psii,15,'linewidth',2); hold on;
[c1,h1]=contour(Z,R,psi_betas4,[-0,0],'--','color',[0 0.4470 0.7410],'linewidth',3);hold on;
contour(Z,R,psi_betas5,[-0,0],'--','color',[0.8500 0.3250 0.0980],'linewidth',3);
contour(Z,R,psi_betas6,[-0,0],'--','color',[0.9290 0.6940 0.1250],'linewidth',3);
contour(Z,R,psi_betas7,[-0,0],'--','color',[0.4940 0.1840 0.5560],'linewidth',3);
plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3); % O-point
set(gca,'Fontsize',fontsize);
text(0.01*zw,R(ido(1)),'O point','fontsize',13);
rm=0.15;zx1=c1(1,2:(end-1));rx1=c1(2,2:(end-1));indx=find(abs(zx1)==min(abs(zx1)));
rs =rx1(indx(1));
xlabel('$z [m]$','interpreter','latex','FontSize',fontsize);
ylabel('$r [m]$','interpreter','latex','FontSize',fontsize);

for jl=1:size(zl,2)
    rectangle('position',[zl(jl)-0.5*wl(jl),rl(jl)-0.5*hl(jl),wl(jl),hl(jl)],...
        'FaceColor', [0.5 0.5 0.5],'linewidth',1);hold on;
end
rectangle('position',[-zw,0,2*zw,rw],'linewidth',1);hold on;
set(gca,'xtick',[-0.75:0.25:0.75]);
text('interpreter','latex','string','(a)','position',[-0.7 0.17],'color','k','FontSize',20)


%%
axes('position',[0.55,0.105,0.425,0.37]);%fig-3
plot(rr*100,Pr_betas4,'-',...
    rr*100,Pr_betas5,rr*100,Pr_betas6,rr*100,Pr_betas7,'linewidth',2);grid on;hold all;
plot([100*rs 100*rs],1.1*[0 1],'--','color','k','LineWidth',1.5);hold all;
set(gca,'Fontsize',fontsize);
xlabel('$r\ [cm], z=0$','interpreter','latex','FontSize',fontsize);
ylabel('$P$','interpreter','latex','FontSize',fontsize);
h=legend('$\beta_s=0.4$','$\beta_s=0.5$','$\beta_s=0.6$','$\beta_s=0.7$');
set(h,'interpreter','latex');legend('boxoff');
xlim([0,100*rm]);set(gca,'xtick',[0:2:100*rm]);set(gca,'ytick',[0:0.2:1.1]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('interpreter','latex','string','$separatrix$','position',[5 0.2],'FontSize',20);
text('interpreter','latex','string','(c)','position',[1 0.95],'color','k','FontSize',20)

axes('position',[0.065,0.105,0.425,0.37]);%fig-2 [left bottom width height]
betas=[0.4, 0.5,0.6,0.7];
n_factor=[1.17,2.12,3.25,3.52];
plot(betas,n_factor,'ko-','linewidth',2,'markersize',8);grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
xlabel('$\beta_s$','interpreter','latex','FontSize',fontsize);
ylabel('$m$','interpreter','latex','FontSize',fontsize);
set(gca,'xtick',[0.4:0.1:0.7]);ylim([0,3.6]);
text('interpreter','latex','string','(b)','position',[0.42 3.05],'color','k','FontSize',20)



