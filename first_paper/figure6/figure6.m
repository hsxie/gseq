close all;clear;clc;
load('MRR2d_betas=0.4_delta=0.148.mat');
bz_4=Bz1d;h_4=h1d;j_4=j1d;p_4=pr1d;psi_4=psi1d;rr_4=rr1d;delta_1=0.148;

load('MRR2d_betas=0.6_delta=0.148.mat');
bz_6=Bz1d;h_6=h1d;j_6=j1d;p_6=pr1d;psi_6=psi1d;rr_6=rr1d;delta_3=0.148;

load('MRR2d_betas=0.8_delta=0.148.mat');
bz_8=Bz1d;h_8=h1d;j_8=j1d;p_8=pr1d;psi_8=psi1d;rr_8=rr1d;delta_5=0.148;

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
plot(rr1d,psi_4,rr1d,psi_6,'k-',...
    rr1d,psi_8,'-','linewidth',3,'markersize',9); hold on;grid on;
set(gca,'Fontsize',fontsize);
ylabel('$\psi$','interpreter','latex','Fontsize',30);
xlabel('$r$','interpreter','latex','Fontsize',35);
text('Interpreter','latex','String','$(c)$','Position',[0.05 0.53],'FontSize',25)
xlim([0,1.5]);set(gca,'ytick',-0.2:0.2:0.6);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
h=legend('$\beta_s=0.4$','$\beta_s=0.6$','$\beta_s=0.8$','location','best');
set(h,'interpreter','latex','fontsize',35);legend('boxoff');

axes(p(1)); 
plot(rr1d,p_4,rr1d,p_6,'k-',rr1d,p_8,'-',...
'linewidth',3,'markersize',9); hold on;grid on;
set(gca,'Fontsize',fontsize);
ylabel('$P$','interpreter','latex','Fontsize',30);grid on;
ylim([0,1.05]);set(gca,'ytick',0.0:0.2:1);set(gca,'xtick',0.0:0.5:2);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$(a)$','Position',[0.05 0.9],'FontSize',25)

axes(p(2)); 
plot(rr1d,j_4,rr1d,j_6,'k-',...
    rr1d,j_8,'-','linewidth',3,'markersize',9); hold on;grid on;
set(gca,'Fontsize',fontsize);
ylabel('$J_\theta$','interpreter','latex','Fontsize',30);grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gca,'ytick',0:5:20);set(gca,'xtick',0.0:0.5:2);
text('Interpreter','latex','String','$(b)$','Position',[0.05 13.7],'FontSize',25)


axes(p(4)); 
plot(rr1d,bz_4,rr1d,bz_6,'k-',...
    rr1d,bz_8,'-','linewidth',3,'markersize',9); hold on;grid on;
set(gca,'Fontsize',fontsize);
ylabel('$B_z$','interpreter','latex','Fontsize',30);grid on;
xlabel('$r$','interpreter','latex','Fontsize',35);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
ylim([-1,1.05]);set(gca,'ytick',-1:0.5:1);set(gca,'xtick',0.0:0.5:2);
text('Interpreter','latex','String','$(d)$','Position',[0.05 0.85],'FontSize',25)