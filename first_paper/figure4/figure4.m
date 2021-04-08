close all;clear;clc;
load('MRR2d_betas=0.7_delta=1.mat');
bz_1=Bz1d;h_1=h1d;j_1=j1d;p_1=pr1d;psi_1=psi1d;rr_1=rr1d;delta_1=1;

load('MRR2d_betas=0.7_delta=0.8.mat');
bz_80=Bz1d;h_80=h1d;j_80=j1d;p_80=pr1d;psi_80=psi1d;rr_80=rr1d;delta_2=0.80;

load('MRR2d_betas=0.7_delta=0.6.mat');
bz_60=Bz1d;h_60=h1d;j_60=j1d;p_60=pr1d;psi_60=psi1d;rr_60=rr1d;delta_3=0.60;

load('MRR2d_betas=0.7_delta=0.4.mat');
bz_40=Bz1d;h_40=h1d;j_40=j1d;p_40=pr1d;psi_40=psi1d;rr_40=rr1d;delta_4=0.40;

load('MRR2d_betas=0.7_delta=0.2.mat');
bz_20=Bz1d;h_20=h1d;j_20=j1d;p_20=pr1d;psi_20=psi1d;rr_20=rr1d;delta_5=0.2;

load('MRR2d_betas=0.7_delta=0.18.mat');
bz_18=Bz1d;h_18=h1d;j_18=j1d;p_18=pr1d;psi_18=psi1d;rr_18=rr1d;delta_6=0.18;


fontsize=30;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.02,0.2,0.5,0.4]);
%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 2;     % 子图行数
TightPlot.RowNumber = 2;    % 子图列数
TightPlot.GapW = 0.06;  % 子图之间的左右间距
TightPlot.GapH = 0.075;   % 子图之间的上下间距
TightPlot.MarginsLower = 0.105;   % 子图与图片下方的间距
TightPlot.MarginsUpper = 0.04;  % 子图与图片上方的间距
TightPlot.MarginsLeft = 0.07;   % 子图与图片左方的间距
TightPlot.MarginsRight = 0.03;  % 子图与图片右方的间距
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]);    % 具体设置参数上一节已经输入了

axes(p(3)); 
plot(rr1d,psi_1,rr1d,psi_60,'k-',...
    rr1d,psi_18,'-','linewidth',3,'markersize',9); hold on;grid on;
set(gca,'Fontsize',fontsize);
ylabel('$\psi$','interpreter','latex','Fontsize',30);
xlabel('$r$','interpreter','latex','Fontsize',35);
text('Interpreter','latex','String','$(c)$','Position',[0.05 0.53],...
'FontSize',25)
xlim([0,1.5]);set(gca,'ytick',-0.2:0.2:0.6);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
h=legend('$\delta_s=1.0\delta_{s}^{RR}$','$\delta_s=0.6\delta_{s}^{RR}$',...
    '$\delta_s=0.18\delta_{s}^{RR}$','location','best');
set(h,'interpreter','latex','fontsize',35);legend('boxoff');

axes(p(1)); 
plot(rr1d,p_1,rr1d,p_60,'k-',rr1d,p_18,'-',...
'linewidth',3,'markersize',9); hold on;grid on;
set(gca,'Fontsize',fontsize);
ylabel('$P$','interpreter','latex','Fontsize',30);
ylim([0,1.05]);set(gca,'ytick',0.0:0.2:1);set(gca,'xtick',0.0:0.5:2);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$(a)$','Position',[0.05 0.9],'FontSize',25)

axes(p(2)); 
plot(rr1d,j_1,rr1d,j_60,'k-',...
    rr1d,j_18,'-','linewidth',3,'markersize',9); hold on;grid on;
set(gca,'Fontsize',fontsize);
ylabel('$J_\theta$','interpreter','latex','Fontsize',30);grid on;
set(gca,'ytick',0:5:20);set(gca,'xtick',0.0:0.5:2);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$(b)$','Position',[0.05 18],'FontSize',25)


axes(p(4)); 
plot(rr1d,bz_1,rr1d,bz_60,'k-',...
    rr1d,bz_18,'-','linewidth',3,'markersize',9); hold on;grid on;
set(gca,'Fontsize',fontsize);
ylabel('$B_z$','interpreter','latex','Fontsize',30);
xlabel('$r$','interpreter','latex','Fontsize',35);
ylim([-0.6,1.05]);set(gca,'ytick',-0.5:0.5:1);set(gca,'xtick',0.0:0.5:2);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$(d)$','Position',[0.05 0.85],'FontSize',25)