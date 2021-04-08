% 19-06-20 14:40
close all;
load('MRR1d.mat');
mu0=4*pi*1e-7;

fontsize=25;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,'defaulttextinterpreter','latex',...
    'position',[0.6,0.2,0.2,0.35],'color', [1, 1, 1]);
% [c1,h1]=contour(Z,R,psi,[-0,0],'r--','linewidth',3);
rx1=c1(2,2:(end-1));zx1=c1(1,2:(end-1));
indx=find(abs(zx1)==min(abs(zx1)));

% axes('position',[0.11,0.1,0.87,0.85]);%fig-1 [left bottom width height]
axes('position',[0.105,0.125,0.85,0.81]);%fig-1 [left bottom width height]
Jzetar=-Jzeta(:,Z0);
Pr=P(:,Z0);
J0=max(abs(Jzetar));
P0=max(abs(Pr));
Bz_mid=Bz(:,Z0);
Bz_2=Bz_mid.^2/2/mu0;
pr_2=Bz_2+Pr;
Be_2=Be^2/2/mu0.*ones(1,size(rr,2));

Br_z=diff(Br,1,2)./(zz(2)-zz(1));
Br_zmid=Br_z(:,Z0);
br_curave=Bz_mid.*Br_zmid/mu0;
sum_br_curave=cumtrapz(rr,br_curave);
plot(rr*100,Pr,'m-','linewidth',3);grid on;hold all;
plot(rr*100,Bz_2,'-','linewidth',3);hold all;
plot(rr*100,pr_2,'r-','linewidth',3);hold all;
plot(rr*100,Be_2,'b--','linewidth',3);hold all;
plot(rr*100,Be_2+sum_br_curave,'--','linewidth',3);hold all;
plot([100*rs 100*rs],[0 5e5],'-.k','LineWidth',1.5);hold all;
set(gca,'Fontsize',fontsize);
xlabel('$r\ (cm), z=0$','interpreter','latex','FontSize',fontsize);
ylabel('$P$','interpreter','latex','FontSize',fontsize);
h=legend('$P$','$B_z^2/(2\mu_0)$','$P+B_z^2/(2\mu_0)$','$B_e^2/(2\\mu_0)$',...
    '$B_e^2/(2\\mu_0)+\int \frac{B_z}{\mu_0}(\frac{\partial{B_r}}{\partial{z}})dr$');
set(h,'interpreter','latex');legend('boxoff');
xlim([0,100*rm]);set(gca,'xtick',[0:2.5:100*rm]);ylim([0,4.5e5]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$separatrix$', 'Position',...
    [100*rs-2 0.4e5],'FontSize',22)





