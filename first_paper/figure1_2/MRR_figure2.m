% 2020-08-02 11:54, Hua-sheng Xie, ENN
% Compare different model, for Tuszewski1984
% datq_mrr0=[rr;Pr;psi;Bz;Jt].';

close all;clear;clc;
load('./RR/rr_datq_t84.mat');
load('./2PE/2pe_datq_t84.mat');
load('./MRR1/mrr1_datq_t84.mat');
load('./MRR2/mrr2_datq_t84.mat');
% load('./MRR2/mrr2_datq.mat');
% load('./3PE/3pe_datq_t84.mat');
load('./3PE/3pe_datq.mat');
load('./dat_t84.mat');

fontsize=30;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.02,0.2,0.5,0.3]);
%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 1;    
TightPlot.RowNumber = 2;   
TightPlot.GapW = 0.05;  
TightPlot.GapH = 0.07;   
TightPlot.MarginsLower = 0.14;   
TightPlot.MarginsUpper = 0.04;  
TightPlot.MarginsLeft = 0.07;   
TightPlot.MarginsRight = 0.02;  
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]);    

axes(p(1));
plot(datq_rr(:,1),datq_rr(:,2),'k-',datq_mrr0(:,1),datq_mrr0(:,2),datq_mrr2(:,1),datq_mrr2(:,2),'-',...
    datq_2pe(:,1),datq_2pe(:,2),'-',datq_3pe(:,1),datq_3pe(:,2),'-',...
    'linewidth',3,'markersize',9); hold on;grid on;
plot(dat_t84(1:3:end,1),dat_t84(1:3:end,2),'r--','markersize',8,'linewidth',2.5);hold on;
set(gca,'Fontsize',fontsize);
ylabel('$P$','interpreter','latex','Fontsize',30);grid on;%xlabel('r');
xlim([0,1.5]);ylim([0,1.05]);
h=legend('$RR$','$MRR-1$','$MRR-2$','$2PE$','$3PE$','$n/n_m$','location','best');
set(h,'interpreter','latex','FontSize',25);legend('boxon');
text('Interpreter','latex','String','$(a)$','Position',[0.05 0.9],'FontSize',25)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gca,'ytick',0.0:0.2:1);set(gca,'xtick',0.0:0.5:1.5);
xlabel('$r$','interpreter','latex','Fontsize',35);

axes(p(2));
plot(datq_rr(:,1),datq_rr(:,5),'k-',datq_mrr0(:,1),datq_mrr0(:,5),datq_mrr2(:,1),datq_mrr2(:,5),'-',...
    datq_2pe(:,1),datq_2pe(:,5),'-',datq_3pe(:,1),datq_3pe(:,5),'-',...
    'linewidth',3,'markersize',9);grid on;
set(gca,'Fontsize',fontsize);
ylabel('$J_\theta$','interpreter','latex','Fontsize',30); hold on;grid on;
xlabel('$r$','interpreter','latex','Fontsize',35);
set(gca,'ytick',0.0:2:10);set(gca,'xtick',0.0:0.5:1.5);
text('Interpreter','latex','String','$(b)$','Position',[0.05 9],'FontSize',25)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);


%% small figure

h2=axes('Position',[0.23 0.6 0.37 0.14],'FontSize',10);
axes(h2);                 
plot(datq_rr(:,1),datq_rr(:,2),'k-',datq_mrr0(:,1),datq_mrr0(:,2),datq_mrr2(:,1),datq_mrr2(:,2),'-',...
    datq_2pe(:,1),datq_2pe(:,2),'-',datq_3pe(:,1),datq_3pe(:,2),'-',...
    'linewidth',3,'markersize',9); hold on;grid on;
plot(dat_t84(1:3:end,1),dat_t84(1:3:end,2),'r--','markersize',8,'linewidth',2.5);hold on;
set(gca,'ytick',0.0:0.25:1);set(gca,'xtick',0.1:0.2:1);
set(h2,'xlim',[0.1 0.301],'FontSize',20);

