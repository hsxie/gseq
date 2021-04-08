% 2020-08-02 11:54, Hua-sheng Xie, ENN
% Compare different model, for Ikeyama09
% datq_mrr0=[rr;Pr;psi;Bz;Jt].';

close all;clear;clc;
load('./RR/rr_datq_i09.mat');
load('./2PE/2pe_datq_i09.mat');
load('./MRR1/mrr1_datq_i09.mat');
load('./MRR2/mrr2_datq_i09.mat');
% load('./MRR2/mrr2_datq.mat');
load('./3PE/3pe_datq_i09.mat');
load('./dat_i09.mat');
dat_i09=diag_data;
dat_i09(:,1)=dat_i09(:,1)/5.61;
dat_i09(:,2)=sqrt(dat_i09(:,2));

%
c=0.58;
delta=0.0051;
b0=1;
a=0.0561;
rr=(0:0.01:1).*2*a;
u=2.*rr.^2/a^2-1;
bi=0.5*c*b0.*(u+u.^3);
be=b0*(1-(1-c).*exp(-(rr-a)/delta));
b=bi.*(u<=1)+be.*(u>1);

p_s_f=1-b.^2/b0^2;
% Jr=-(diff(p_s_f))./diff(psi).*(rr(2:end));Jr=[Jr(1);Jr];


fontsize=30;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.02,0.2,0.5,0.3]);

%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 1;     % 子图行数
TightPlot.RowNumber = 2;    % 子图列数
TightPlot.GapW = 0.05;  % 子图之间的左右间距
TightPlot.GapH = 0.07;   % 子图之间的上下间距
TightPlot.MarginsLower = 0.14;   % 子图与图片下方的间距
TightPlot.MarginsUpper = 0.04;  % 子图与图片上方的间距
TightPlot.MarginsLeft = 0.07;   % 子图与图片左方的间距
TightPlot.MarginsRight = 0.02;  % 子图与图片右方的间距
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]);    % 具体设置参数上一节已经输入了

axes(p(1)); 
plot(datq_rr(:,1),datq_rr(:,2).^2,'k-',datq_mrr0(:,1),datq_mrr0(:,2).^2,datq_mrr2(:,1),datq_mrr2(:,2).^2,'-',...
    datq_2pe(:,1),datq_2pe(:,2).^2,'-',datq_3pe(:,1),datq_3pe(:,2).^2,'-',...
    rr/a,p_s_f.^2,'-',...
    'linewidth',3,'markersize',9); hold on;grid on;
plot(dat_i09(1:1:end,1),dat_i09(1:1:end,2).^2,'ro-','markersize',8,'linewidth',2.5);hold on;
set(gca,'Fontsize',fontsize);
ylabel('$P^2$','interpreter','latex','Fontsize',30);grid on;
xlim([0,1.5]);ylim([0,1.05]);
hc=legend('$\mathrm{RR}$','$\mathrm{MRR-1}$','$\mathrm{MRR-2}$','$\mathrm{2PE}$',...
    '$\mathrm{3PE}$','$\mathrm{s-f}$','$\mathrm{i/i_0}$','location','best');
set(hc,'interpreter','latex','FontSize',25);legend('boxon');
text('Interpreter','latex','String','$(a)$','Position',[0.05,0.9],'FontSize',25)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gca,'ytick',0.0:0.2:1);set(gca,'xtick',0.0:0.5:1.5);
xlabel('$r$','interpreter','latex','Fontsize',35);


axes(p(2)); 
plot(datq_rr(:,1),datq_rr(:,5),'k-',datq_mrr0(:,1),datq_mrr0(:,5),datq_mrr2(:,1),datq_mrr2(:,5),'-',...
     datq_2pe(:,1),datq_2pe(:,5),'-',datq_3pe(:,1),datq_3pe(:,5),'-',...
   'linewidth',3,'markersize',9);
set(gca,'Fontsize',fontsize);
ylabel('$J_\theta$','interpreter','latex','Fontsize',30); hold on;grid on;
xlabel('$r$','interpreter','latex','Fontsize',35);
% text('Interpreter','latex','String',0.05,7.0,'(c)','Fontsize',20);
text('Interpreter','latex','String','$(b)$','Position',[0.05,7.0],'FontSize',25)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gca,'ytick',0.0:2:8);set(gca,'xtick',0.0:0.5:1.5);


%% small figure
% 在原图上插入一个新的小图像

h2=axes('Position',[0.23 0.6 0.37 0.14],'FontSize',10);
axes(h2);                 
plot(datq_rr(:,1),datq_rr(:,2).^2,'k-',datq_mrr0(:,1),datq_mrr0(:,2).^2,datq_mrr2(:,1),datq_mrr2(:,2).^2,'-',...
    datq_2pe(:,1),datq_2pe(:,2).^2,'-',datq_3pe(:,1),datq_3pe(:,2).^2,'-',...
    rr/a,p_s_f.^2,'-','linewidth',3);hold on;
plot(dat_i09(1:1:end,1),dat_i09(1:1:end,2).^2,'ro-','markersize',8,'linewidth',2.5);hold on;
set(gca,'ytick',0.0:0.25:1);set(gca,'xtick',0.0:0.1:1);
set(h2,'xlim',[0.3 0.4],'FontSize',20);



