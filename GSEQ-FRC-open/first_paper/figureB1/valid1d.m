% 2020-08-04 07:46, Hua-sheng Xie, ENN
% Valid the 1D equilibrium use Hill vortex psi;

close all;clear;clc;

Be=1.0;
a=1; % a=rs
dat=[1,0.1;
    1,0.5;
    1,0.9;
    2,0.1;
    2,0.9;
    3,0.5;
    3,0.9;
    4,0.5;];


fontsize=30;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.02,0.2,0.25,0.3]);
%% INPUT TIGHTPLOT PARAMETERS
TightPlot.ColumeNumber = 1;     % 子图行数
TightPlot.RowNumber = 1;    % 子图列数
TightPlot.GapW = 0.09;  % 子图之间的左右间距
TightPlot.GapH = 0.05;   % 子图之间的上下间距
TightPlot.MarginsLower = 0.15;   % 子图与图片下方的间距
TightPlot.MarginsUpper = 0.04;  % 子图与图片上方的间距
TightPlot.MarginsLeft = 0.15;   % 子图与图片左方的间距
TightPlot.MarginsRight = 0.03;  % 子图与图片右方的间距
p = tight_subplot(TightPlot.ColumeNumber,TightPlot.RowNumber,...
    [TightPlot.GapH TightPlot.GapW],...
    [TightPlot.MarginsLower TightPlot.MarginsUpper],...
    [TightPlot.MarginsLeft TightPlot.MarginsRight]);    % 具体设置参数上一节已经输入了

nj=size(dat,1);
for jd=1:nj
%     E=1; % elongation
%     betas=0.3;
    
    E=dat(jd,1);
    betas=dat(jd,2);
    
    b=a*E;
    z=0;
    
    dr=0.01;
    rr=0:dr:1.0;
    
    B0=sqrt(1-betas);
    psi=B0/2*rr.^2.*(1-rr.^2/a^2-z^2/b^2);
    Br=B0*rr*z/b^2;
    Bz=B0.*(1-2*rr.^2/a^2-z^2/b^2);
    dBrdz=B0*rr/b^2;
    
    Pm=1+cumsum(2*Bz.*dBrdz*dr);
    
    str{jd}=['$E=$',num2str(E),'$, \beta_s=$',num2str(betas)];
    
    plot(rr,(Pm-Pm(1))/(Pm(1)),'linewidth',3); hold on;grid on;
%     ylim([0,1.3]);
end
legend(str,'location','best','interpreter','latex');
legend('boxoff');
ylabel('$(P_m-P_{m0})/P_{m0}$','interpreter','latex','Fontsize',25);
xlabel('$r$','interpreter','latex','Fontsize',35);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gca,'xtick',[0:0.2:1]);


% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
%   'PaperSize',[screenposition(3:4)]);
% print(gcf,'-dpng',['valid_1d.png']);
