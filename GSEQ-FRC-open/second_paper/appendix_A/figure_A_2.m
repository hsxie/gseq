% 19-06-20 14:40
close all;
fontsize=25;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,'defaulttextinterpreter','latex',...
    'position',[0.5,0.2,0.28,0.6],'color', [1, 1, 1]);
load('MRR1d.mat');
mu0=4*pi*1e-7;

% axes('position',[0.11,0.75,0.85,0.23]);%fig-1
% axes('position',[0.065,0.69,0.91,0.3]);%fig-1[left bottom width height]
axes('position',[0.11,0.75,0.85,0.23]);%fig-1
ido=find(psi==min(min(psi)));
indtmp= psi>0; psio=psi; psio(indtmp)=NaN;
indtmp=find(psi<0); psii=psi; psii(indtmp)=NaN;
contour(Z,R,psio,7,'linewidth',2); hold on;
contour(Z,R,psii,30,'linewidth',2); hold on;%colormap('hot'); 
[c1,h1]=contour(Z,R,psi,[-0,0],'r--','linewidth',3);
plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3); % O-point
text(0.01*zw,R(ido(1)),'O point','fontsize',13);

ylabel('r [m]');
zx1=c1(1,2:(end-1)); rx1=c1(2,2:(end-1)); zx_no=size(zx1,2);
indx=find(abs(zx1)==min(abs(zx1)));
set(gca,'ytick',[0:0.05:0.15]);
set(gca,'xtick',[-0.6:0.3:0.6]);
min(psi(:,Z0));
rs =rx1(indx(1));
plot([-0.75,0.75],[0.14,0.14],'k-',[-0.75,-0.75],[0,0.14],'k-',[0.75,0.75],[0,0.14],...
    'k-','linewidth',3);hold on;
h = colorbar('east');
% h.TickLabels = (0:1:6);
% h.Label.String = '\times 10^{-2}';

h=rectangle('Position',[-0.75 0.147 0.1 0.006],'FaceColor', [0.5 0.5 0.5],'linewidth',1);
w = h.LineWidth;h.LineWidth = 3;hold on;
h=rectangle('Position',[-0.60 0.147 1.2 0.006],'FaceColor', [0.5 0.5 0.5],'linewidth',1);
w = h.LineWidth;h.LineWidth = 3;hold on;
h=rectangle('Position',[0.65 0.147 0.1 0.006],'FaceColor', [0.5 0.5 0.5],'linewidth',1);
w = h.LineWidth;h.LineWidth = 3;hold on;

ylabel('$r\ [m]$','interpreter','latex','FontSize',fontsize);
xlabel('$z\ [m]$','interpreter','latex','FontSize',fontsize);
text('Interpreter','latex','String','$(a)$','Position',[-0.7 0.12],'color','k',...
    'FontSize',fontsize)

% axes('position',[0.12,0.41,0.38,0.26]); %fig-2 [left bottom width height]
% axes('position',[0.065,0.39,0.425,0.24]);%fig-2[left bottom width height]
axes('position',[0.12,0.41,0.38,0.26]); %fig-2 [left bottom width height]
psiw_coil=psiw_coil.*CI;
psiw_coil_u=psiw_coil(nR:nR+nZ-1);
plot(Z(end,:),psiw_coil_u,'r-',Z(end,:),psi(end,:),'g--','linewidth',3);grid on;
% plot(zb,psib,'r-',Z(1,:),psi_plasma,'g--','linewidth',3);grid on;hold all;
set(gca,'Fontsize',fontsize);
xlabel('$z [m]$','interpreter','latex','FontSize',fontsize);
ylabel('$\psi \ at \ wall$','interpreter','latex','FontSize',fontsize);
xlim([-0.75,0.75]);set(gca,'xtick',[-0.6:0.4:0.6]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gcf, 'PaperPositionMode', 'auto');
text('Interpreter','latex','String','$(b)$','Position',[-0.7 0.04],'color','k',...
    'FontSize',fontsize);
h=legend('$\mathbf{w/o \ plasma}$','$\mathbf{w \ plasma}$');set(h,'interpreter','latex');
legend('boxoff');

% axes('position',[0.58,0.41,0.38,0.26]); %fig-3
% axes('position',[0.55,0.39,0.425,0.24]);%fig-3[left bottom width height]
axes('position',[0.58,0.41,0.38,0.26]); %fig-3
Bz_coil_u=Bz_coil.*CI;
Bz_coil_u=Bz_coil_u(nR:nR+nZ-1);
plot(Z(end,:),Bz_coil_u,'r-',Z(end,:),Bz(end,:),'g--','linewidth',3);grid on;hold all;
ylabel('$B_{z} \ at \ wall$','interpreter','latex','FontSize',fontsize);
xlabel('$z [m]$','interpreter','latex','FontSize',fontsize);
xlim([-0.75,0.75]);set(gca,'xtick',[-0.6:0.4:0.6]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$(c)$','Position',[-0.7 5.6],'color','k','FontSize',fontsize)
h=legend('$\mathbf{w/o \ plasma}$','$\mathbf{w \ plasma}$');set(h,'interpreter','latex');
legend('boxoff');


% axes('position',[0.12,0.07,0.38,0.26]); %fig-4 [left bottom width height]
% axes('position',[0.065,0.065,0.425,0.24]);%fig-4[left bottom width height]
axes('position',[0.12,0.07,0.38,0.26]); %fig-4 [left bottom width height]
Jzetar=-Jzeta(:,Z0);
Pr=P(:,Z0);
J0=max(abs(Jzetar));
P0=max(abs(Pr));
plot(rr*100,Pr,'m-',rr1d.*Rs*100,pr1d.*(Be^2/(2*mu0)),'--','linewidth',3);grid on;hold all;
plot([100*rs 100*rs],[0 5e5],'-.k','LineWidth',1.5);hold all;
set(gca,'Fontsize',fontsize);
xlabel('$r\ [cm], z=0$','interpreter','latex','FontSize',fontsize);
ylabel('$P$','interpreter','latex','FontSize',fontsize);
h=legend('$2D$','$1D$');set(h,'interpreter','latex');legend('boxoff');
xlim([0,100*rm]);set(gca,'xtick',[0:2.5:100*rm]);
ylim([0,4.3e5]);set(gca,'ytick',[0:1e5:4.5e5]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$(d)$','Position',[1 3.8e5],'color','k','FontSize',fontsize);
text('Interpreter','latex','String','$separatrix$', 'Position',...
    [100*rs-2 1e5],'FontSize',22)


% axes('position',[0.58,0.07,0.38,0.26]); %fig-5
% axes('position',[0.55,0.065,0.425,0.24]);%fig-5[left bottom width height]
axes('position',[0.58,0.07,0.38,0.26]); %fig-5
Jzetar=-Jzeta(:,Z0);
plot(rr*100,Jzetar,'m-',rr1d.*Rs*100,j1d.*(Be/(2*mu0*Rs)),'--','linewidth',3);grid on;hold all;
plot([100*rs 100*rs],1.1*[0 J0],'-.k','LineWidth',1.5);hold all;
set(gca,'Fontsize',fontsize);
xlabel('$r\ [cm], z=0$','interpreter','latex','FontSize',fontsize);
ylabel('$J_{\theta}$','interpreter','latex','FontSize',fontsize);
h=legend('$2D$','$1D$');set(h,'interpreter','latex');legend('boxoff');
xlim([0,100*rm]);set(gca,'xtick',[0:2.5:100*rm]);
ylim([0,5e7]);set(gca,'ytick',[0:1e7:5e7]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gcf, 'PaperPositionMode', 'auto');
text('Interpreter','latex','String','$(e)$','Position',[1 4.5e7],'color','k','FontSize',fontsize);


% figure(2)
% plot(rr*100,Bz(:,floor(nZ/2)))

