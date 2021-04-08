% 19-06-20 14:40
close all;clc;
fontsize=25;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,'defaulttextinterpreter','latex',...
    'position',[0.5,0.2,0.28,0.6],'color', [1, 1, 1]);%0.02,0.2,0.28,0.5
Bz_plasma=Bz(end,:);
psi_plasma=psi(end,:);

load('MRR1d_beta=0_6_delta=0_4.mat');
mu0=4*pi*1e-7;

axes('position',[0.11,0.75,0.85,0.23]);%fig-1
ido=find(psi==min(min(psi)));
indtmp= psi>0; psio=psi; psio(indtmp)=NaN;
indtmp=find(psi<0); psii=psi; psii(indtmp)=NaN;
contour(Z,R,psio,5,'linewidth',2); hold on;
contour(Z,R,psii,30,'linewidth',2); hold on;colorbar('east');
[c1,h1]=contour(Z,R,psi,[-0,0],'r--','linewidth',3);
plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3); hold on;% O-point
text(0.01*zw,R(ido(1)),'O point','fontsize',13);
plot([min(c1(1,:)),max(c1(1,:))],[0,0],'r--','linewidth',4);hold on;

zx1=c1(1,2:(end-1)); rx1=c1(2,2:(end-1)); zx_no=size(zx1,2);
indx=find(abs(zx1)==min(abs(zx1)));
strp='p(\psi)=\beta_sexp(\alpha\psi)*(\sigma\psi+1)^2';
xlabel('$z [m]$','interpreter','latex','FontSize',25);
ylabel('$r [m]$','interpreter','latex','FontSize',fontsize);
set(gca,'ytick',[0:0.05:0.15]);set(gca,'xtick',[-0.75:0.25:0.75]);
min(psi(:,Z0));
rs =rx1(indx(1));
plot(zzb,rrb,'k-','linewidth',3);

h=rectangle('Position',[-0.75 0.147 0.1 0.006],'FaceColor', [0.5 0.5 0.5],'linewidth',1);w = h.LineWidth;h.LineWidth = 3;hold on;
h=rectangle('Position',[-0.60 0.147 1.2 0.006],'FaceColor', [0.5 0.5 0.5],'linewidth',1);w = h.LineWidth;h.LineWidth = 3;hold on;
h=rectangle('Position',[0.65 0.147 0.1 0.006],'FaceColor', [0.5 0.5 0.5],'linewidth',1);w = h.LineWidth;h.LineWidth = 3;hold on;
text('Interpreter','latex','String','$(a)$','Position',[-0.7 0.12],'color','k','FontSize',fontsize)

axes('position',[0.12,0.07,0.38,0.26]); %fig-4 [left bottom width height]
Jzetar=-Jzeta(:,Z0);Pr=P(:,Z0);
plot(rr*100,Pr,'m-',rr1d.*Rs*100,pr1d.*(Be^2/(2*mu0)),'--','linewidth',3);grid on;hold all;
plot([100*rs 100*rs],[0 5e5],'-.k','LineWidth',1.5);hold all;
set(gca,'Fontsize',fontsize);
xlabel('$r\ [cm], z=0$','interpreter','latex','FontSize',fontsize);
ylabel('$P$','interpreter','latex','FontSize',fontsize);
h=legend('$2D$','$1D$');set(h,'interpreter','latex');legend('boxoff');
xlim([0,100*rm]);set(gca,'xtick',[0:2.5:100*rm]);
ylim([0,4.5e5]);set(gca,'ytick',[0:1e5:4e5]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$Separatrix$', 'Position',...
    [100*rs-2 1e5],'FontSize',20)
text('Interpreter','latex','String','$(d)$','Position',[1 3.9e5],'color','k','FontSize',fontsize)

axes('position',[0.58,0.07,0.38,0.26]); %fig-5
Jzetar=-Jzeta(:,Z0);J0=max(-Jzeta(:,Z0));
plot(rr*100,-Jzetar,'m',rr1d.*Rs*100,j1d.*(Be/(2*mu0*Rs)),'--','linewidth',3);grid on;hold all;
plot([100*rs 100*rs],[0 5e7],'-.k','LineWidth',1.5);hold all;
set(gca,'Fontsize',fontsize);
xlabel('$r\ [cm], z=0$','interpreter','latex','FontSize',fontsize);
ylabel('$J_{\theta}$','interpreter','latex','FontSize',fontsize);
h=legend('$2D$','$1D$');set(h,'interpreter','latex');legend('boxoff');
xlim([0,100*rm]);set(gca,'xtick',[0:2.5:100*rm]);
ylim([0,5.1e7]);set(gca,'ytick',[0:1e7:4e7]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gcf, 'PaperPositionMode', 'auto');
text('Interpreter','latex','String','$(e)$','Position',[1 4.5e7],'color','k','FontSize',fontsize)

load('psib_paper_beta=0_6_delta=0_4.mat');% read data
axes('position',[0.12,0.41,0.38,0.26]); %fig-2 [left bottom width height]
plot(zb,psib,'r-',Z(1,:),psi_plasma,'g--','linewidth',3);grid on;hold all;
set(gca,'Fontsize',fontsize);
xlabel('$z [m]$','interpreter','latex','FontSize',fontsize);
ylabel('$\psi \ at \ wall$','interpreter','latex','FontSize',fontsize);
xlim([-0.75,0.75]);set(gca,'xtick',[-0.6:0.4:0.6]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gcf, 'PaperPositionMode', 'auto');
h=legend('$\mathbf{w/o \ plasma}$','$\mathbf{w \ plasma}$');set(h,'interpreter','latex');
legend('boxoff');
text('Interpreter','latex','String','$(b)$','Position',[-0.7 0.035],'color','k','FontSize',fontsize);
ylim([0,0.04]);set(gca,'ytick',[0:0.01:0.04]);

axes('position',[0.58,0.41,0.38,0.26]); %fig-3
plot(zb,Bz,'r-',Z(1,:),Bz_plasma,'g--','linewidth',3);grid on;hold all;
ylabel('$B_{z} \ at \ wall$','interpreter','latex','FontSize',fontsize);
xlabel('$z [m]$','interpreter','latex','FontSize',fontsize);
xlim([-0.75,0.75]);set(gca,'xtick',[-0.6:0.4:0.6]);set(gca,'ytick',[0:2:10]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
h=legend('$\mathbf{w/o \ plasma}$','$\mathbf{w \ plasma}$');set(h,'interpreter','latex');
legend('boxoff');
text('Interpreter','latex','String','$(c)$','Position',[-0.7 9],'color','k','FontSize',fontsize)


