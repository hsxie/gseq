close all;clear;clc;
fname='data2.dat';
fid=fopen(fname,'r');

for i=1:210
    betas(i)=fscanf(fid,'%f',1);
    ls_1(i)=fscanf(fid,'%f',1);
    ls_0(i)=fscanf(fid,'%f',1);
    alpha(i)=fscanf(fid,'%f',1);
    sigma(i)=fscanf(fid,'%f',1);
    hx(i)=fscanf(fid,'%f',1);
    psim(i)=fscanf(fid,'%f',1);
    ave_beta(i)=fscanf(fid,'%f',1);
    tot_I(i)=fscanf(fid,'%f',1);
    tot_Ix(i)=fscanf(fid,'%f',1);
    p(i)=fscanf(fid,'%f',1);
    psig(i)=fscanf(fid,'%f',1);
    x0(i)=fscanf(fid,'%f',1);
end
ls=ls_1.*ls_0;

[xq,yq] = meshgrid(min(betas):0.0002*(max(betas)-min(betas)):max(betas),...
    min(ls):0.0002*(max(ls)-min(ls)):max(ls));
% [xq,yq] = meshgrid(min(betas):0.01*(max(betas)-min(betas)):max(betas),...
%     min(ls):0.01*(max(ls)-min(ls)):max(ls));

psim_new = griddata(betas,ls,psim,xq,yq);
alpha_new = griddata(betas,ls,alpha,xq,yq);
sigma_new = griddata(betas,ls,sigma,xq,yq);
hx_new = griddata(betas,ls,hx,xq,yq);
ave_beta_new = griddata(betas,ls,ave_beta,xq,yq);
p_new = griddata(betas,ls,p,xq,yq);

tot_I_new = griddata(betas,ls,tot_I,xq,yq);
tot_I_total = griddata(betas,ls,tot_Ix,xq,yq);
tot_I_1 = griddata(betas,ls,tot_I./tot_Ix,xq,yq);

for i=1:size(yq,1)
    alpha0(i)=8*log((sqrt(1-xq(1,i))+1)/sqrt(xq(1,i)));
    ls0(i)=1./(alpha0(i).*sqrt(1-xq(1,i)));
end

for i=1:size(yq,1)
    for j=1:size(yq,1)
       if yq(j,i) >ls0(i)
           psim_new(j,i)=nan;
           alpha_new(j,i)=nan;
           sigma_new(j,i)=nan;
           hx_new(j,i)=nan;
           ave_beta_new(j,i)=nan;
           tot_I_new(j,i)=nan;
           tot_I_total(j,i)=nan;
           tot_I_1(j,i)=nan;
           p_new(j,i)=nan;
       end
    end
end

fontsize=25;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,'defaulttextinterpreter','latex',...
    'position',[0.02,0.2,0.45,0.38],'color', [1, 1, 1]); %0.02,0.2,0.4,0.38 %0.02,0.2,0.4,0.3
set(gcf, 'PaperPositionMode', 'auto');% 使print出来的与屏幕显示大小相同

bb=0.05:0.02:0.95;
ll=0.*bb;
aa=0.*bb;
pp=0.*bb;
for jb=1:length(bb)
    betas=bb(jb);
    alpha0=8*log((sqrt(1-betas)+1)/sqrt(betas));
    ls0=1/(alpha0*sqrt(1-betas));
    psim0=log(betas)/alpha0;
    aa(jb)=alpha0;
    pp(jb)=psim0;
    ll(jb)=ls0;
end

No=30;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
up=0.565; down=0.13; width1=0.225; hight1=0.42;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes('position',[0.06,up,width1,hight1]);%fig-1 [left bottom width height]
[c,h]=contourf(xq,yq,psim_new,No);hold on;
set(gca,'Fontsize',fontsize);
set(h.Parent,'YScale','log');
ax.XAxis.Visible = 'off';
ylabel('$\delta_s$','fontsize',28,'position',[-0.08 0.35]);
semilogy(bb,ll,'linewidth',2);
ylim([0,max(max(yq))]);set(gca,'Ytick',[0 0.01 0.1  1]);
set(gca,'xticklabel',get(gca,'xtick'),'yticklabel',get(gca,'ytick'));
set(gca,'xtick',[]);
text('Interpreter','latex','String','$(a) \quad \psi_{min}$','Position',[0.07 1.7],'FontSize',25);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);

c = colorbar('Location','east');
% c = colorbar('position',[0.27 0.45 0.02 0.15]);
ax = gca;
axpos = ax.Position;
c.Position(3) = 0.3*c.Position(3);
HorizontalAlignment='right';
ax.Position = axpos;


axes('position',[0.297,up,width1,hight1]);%fig-2 [left bottom width height]
[c,h]=contourf(xq,yq,alpha_new,No);hold on;
semilogy(bb,ll,'linewidth',2);
set(gca,'Fontsize',fontsize);
set(h.Parent,'YScale','log');
ylim([0,max(max(yq))]);set(gca,'Ytick',[0 0.1 1]);
set(gca,'xtick',[],'ytick',[]);
text('Interpreter','latex','String','$(b) \quad \alpha$','Position',[0.07 1.7],'FontSize',25);

c = colorbar('Location','east');
ax = gca;
axpos = ax.Position;
c.Position(3) = 0.3*c.Position(3);
ax.Position = axpos;

axes('position',[0.535,up,width1,hight1]);%fig-3 [left bottom width height]
[c,h]=contourf(xq,yq,sigma_new,No);hold on;
semilogy(bb,ll,'linewidth',2);
set(gca,'Fontsize',fontsize);
set(h.Parent,'YScale','log');
set(gca,'xtick',[],'ytick',[]);
text('Interpreter','latex','String','$(c) \quad \sigma$','Position',[0.07 1.7],'FontSize',25);

c = colorbar('Location','east','Ticks',[2,8,14]);
ax = gca;
axpos = ax.Position;
c.Position(3) = 0.3*c.Position(3);
ax.Position = axpos;

axes('position',[0.77,up,width1,hight1]);%fig-4 [left bottom width height]
[c,h]=contourf(xq,yq,hx_new,No);hold on;
semilogy(bb,ll,'linewidth',2);
set(gca,'Fontsize',fontsize);
set(h.Parent,'YScale','log');
xlabel('$\beta_s$');
set(gca,'xtick',[],'ytick',[]);
text('Interpreter','latex','String','$(d) \quad h$','Position',[0.07 1.7],'FontSize',25);

c = colorbar('Location','east','Ticks',[0.5,1,1.5]);
ax = gca;
axpos = ax.Position;
c.Position(3) = 0.3*c.Position(3);
ax.Position = axpos;

axes('position',[0.06,down,width1,hight1]);%fig-5 [left bottom width height]
[c,h]=contourf(xq,yq,ave_beta_new,No);hold on;
semilogy(bb,ll,'linewidth',2);
set(gca,'Fontsize',fontsize);
set(h.Parent,'YScale','log');
xlabel('$\beta_s$','fontsize',28);ylabel('$\delta_s$','fontsize',28,'position',[-0.08 0.35]);
ylim([0,max(max(yq))]);set(gca,'Ytick',[0 0.1 1]);set(gca,'Xtick',[0.2 0.5 0.8]);
set(gca,'xticklabel',get(gca,'xtick'),'yticklabel',get(gca,'ytick'));
text('Interpreter','latex','String','$(e) \quad \left \langle \beta \right \rangle$',...
    'Position',[0.07 1.7],'FontSize',20);

c = colorbar('Location','east','Ticks',[0.5,0.7,0.9]);
ax = gca;
axpos = ax.Position;
c.Position(3) = 0.3*c.Position(3);
ax.Position = axpos;

axes('position',[0.297,down,width1,hight1]);%fig-6 [left bottom width height]
[c,h]=contourf(xq,yq,tot_I_new,No);hold on;
semilogy(bb,ll,'linewidth',2);
set(gca,'Fontsize',fontsize);
set(h.Parent,'YScale','log');
xlabel('$\beta_s$','fontsize',28);set(gca,'ytick',[]);set(gca,'Xtick',[0.2 0.5 0.8]);
text('Interpreter','latex','String','$(f) \quad I_{in}$','Position',[0.07 1.7],'FontSize',25);

c = colorbar('Location','east');
ax = gca;
axpos = ax.Position;
c.Position(3) = 0.3*c.Position(3);
ax.Position = axpos;

axes('position',[0.535,down,width1,hight1]);%fig-7 [left bottom width height]
[c,h]=contourf(xq,yq,tot_I_1,No);hold on;
semilogy(bb,ll,'linewidth',2);
set(gca,'Fontsize',fontsize);
set(h.Parent,'YScale','log');
xlabel('$\beta_s$','fontsize',28);set(gca,'ytick',[]);set(gca,'Xtick',[0.2 0.5 0.8]);
text('Interpreter','latex','String','$(g) \quad c_{i}$','Position',[0.07 1.7],'FontSize',25);

c = colorbar('Location','east','Ticks',[0.5,0.7,0.9]);
ax = gca;
axpos = ax.Position;
c.Position(3) = 0.3*c.Position(3);
ax.Position = axpos;

axes('position',[0.77,down,width1,hight1]);%fig-8 [left bottom width height]
[c,h]=contourf(xq,yq,p_new,No);hold on;
semilogy(bb,ll,'linewidth',2);
set(gca,'Fontsize',fontsize);
set(h.Parent,'YScale','log');
xlabel('$\beta_s$','fontsize',28);set(gca,'ytick',[]);set(gca,'Xtick',[0.2 0.5 0.8]);
text('Interpreter','latex','String','$(h) \quad c_{p}$','Position',[0.07 1.7],'FontSize',25);

c = colorbar('Location','east','Ticks',[0.5,0.7,0.9]);
ax = gca;
axpos = ax.Position;
c.Position(3) = 0.3*c.Position(3);
ax.Position = axpos;
colormap('hot');