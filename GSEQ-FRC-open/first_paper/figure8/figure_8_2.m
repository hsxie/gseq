close all;
fontsize=25;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,'defaulttextinterpreter','latex',...
    'position',[0.5,0.4,0.5,0.4],'color', [1, 1, 1]);%0.02,0.2,0.28,0.5
%0.55,0.2,0.4,0.5
set(gcf, 'PaperPositionMode', 'auto');% 使print出来的与屏幕显示大小相同

C1=['m','k','b','g','c','r'];
i=1;
load('gs1d.mat');

% subplot(3,2,1:2);
axes('position',[0.08,0.74,0.91,0.25]);%fig-1[left bottom width height]
ido=find(psi==min(min(psi)));
indtmp= psi>0; psio=psi; psio(indtmp)=NaN;
indtmp=find(psi<0); psii=psi; psii(indtmp)=NaN;
contour(Z,R,psio,8,'linewidth',2); hold on;
contour(Z,R,psii,15,'linewidth',2); hold on;%colorbar('northOutside');
[c1,h1]=contour(Z,R,psi,[-0,0],'r--','linewidth',3);
plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3); % O-point
set(gca,'Fontsize',fontsize);
text(0.01*zw,R(ido(1)),'O point','fontsize',13);

ylabel('$r (m)$','Fontsize',fontsize);
zx1=c1(1,2:(end-1)); rx1=c1(2,2:(end-1)); zx_no=size(zx1,2);
indx=find(abs(zx1)==min(abs(zx1)));
ylim([0,0.18]);
set(gca,'ytick',[0:0.05:0.22]);
% min(psi(:,floor(Zs*nZ)))
rs =rx1(indx(1));
plot(zzb,rrb,'k-','linewidth',3);
xlabel('$Z(m)$','position',[-0.2 -0.025],'Fontsize',fontsize);
text('Interpreter','latex','String','$(a)$','Position',[-0.7 0.163 1],'FontSize',fontsize)

Z_mid=floor(nZ/2)+1;
Pr=P(:,Z_mid);
Bzr=Bz(:,Z_mid);
Jzetar=Jzeta(:,Z_mid);
psir=psi(:,Z_mid);

% subplot(3,2,3);
axes('position',[0.08,0.41,0.425,0.25]);%fig-2[left bottom width height]
plot(100.*r_rr,pmrr/max(pmrr),'-',100.*rr,Pr/max(Pr),'--',100.*r_rr,prr/max(prr),'LineWidth',3);hold on;
set(gca,'Fontsize',fontsize);grid on;
ylabel('$P/P_{m}$','FontSize',fontsize);
h=legend('$MRR-1D$','$MRR-2D$','$RR$','location','best');set(h,'interpreter','latex');
legend('boxoff');
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gca,'Ytick',[0:0.25:1]);
text('Interpreter','latex','String','$(b)$','Position',[0.2 0.8],'FontSize',fontsize)
xlim([0,100*rm])

% subplot(3,2,4);
axes('position',[0.565,0.41,0.425,0.25]);%fig-3[left bottom width height]
plot(100.*r_rr,Bzmrr,'-',100.*rr,Bzr,'--',100.*r_rr,Bzrr,'LineWidth',3);hold on;
set(gca,'Fontsize',25);grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
ylabel('$B_z (T)$','Fontsize',fontsize);
text('Interpreter','latex','String','$(c)$','Position',[0.2 0.35],'FontSize',fontsize)
xlim([0,100*rm]);set(gca,'Xtick',[0:5:20]);

% subplot(3,2,5);
axes('position',[0.08,0.095,0.425,0.25]);%fig-4[left bottom width height]
plot(100.*r_rr,Jr_mrr/max(Jr_mrr),'-',...
    100.*rr,-Jzetar/max(-Jzetar),'--',100.*r_rr,Jr_rr/max(Jr_rr),'LineWidth',3);hold on;
set(gca,'Fontsize',fontsize);grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
ylabel('$J_{\theta}/ J_m$','Fontsize',fontsize);set(gca,'Ytick',[0:0.25:1]);
text('Interpreter','latex','String','$(d)$','Position',[0.2 0.8],'FontSize',fontsize)
xlim([0,100*rm]);
xlabel('$r [cm] \ with \ Z=0$','interpreter','latex','Fontsize',fontsize);


% subplot(3,2,6);
axes('position',[0.565,0.095,0.425,0.25]);%fig-5[left bottom width height]
rr_data=importdata('rr_data.mat',' ');
rw1=rr_data(:,1);P2w1=rr_data(:,2);
diag_data=importdata('diag_data.mat',' ');
rw2=diag_data(:,1);P2w2=diag_data(:,2);
f_data=importdata('f.mat',' ');
rf=f_data(:,1);P2f=f_data(:,2);
sf_data=importdata('s_f.mat',' ');
rsf=sf_data(:,1);P2sf=sf_data(:,2);
p_data=importdata('p.mat',' ');
rp=p_data(:,1);P2p=p_data(:,2);
load('mrr.mat');
mrr1d_r=rx;mrr1d_p2=Px2;
plot(100.*r_rr,psi_r/abs(min(psi_r)),100.*rr,psir/abs(min(psir)),'--',...
    'LineWidth',3);hold on;
set(gca,'Fontsize',fontsize);grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
xlim([0,6]);set(gca,'Ytick',[-1:0.5:1]);
xlabel('$r [cm] \ with \ Z=0$','interpreter','latex','Fontsize',fontsize);
ylabel('$\psi / \psi_{min}$','interpreter','latex','Fontsize',fontsize);
text('Interpreter','latex','String','$(e)$','Position',[0.08 0.6],'FontSize',fontsize);
