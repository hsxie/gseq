close all;clear;clc;
fontsize=25;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,'defaulttextinterpreter','latex',...
    'position',[0.5,0.2,0.28,0.6],'color', [1, 1, 1]);%0.02,0.2,0.28,0.5

load('psi.mat');
psi_matlab =psiq;

load('psi_grass1.mat');
psi_grass =Z1;

psi =psi_grass-psi_matlab;

rmin =0;    zmin =-1.5;     rmax =0.17;     zmax =1.5;

r =linspace(rmin,rmax,size(psi,1));
z =linspace(zmin,zmax,size(psi,2));
[X,Y]=ndgrid(r,z);

rr =linspace(rmin,rmax,size(psi,1));
zz =linspace(zmin,zmax,size(psi,2));
[X1,Y1] =ndgrid(rr,zz);

load('psib.mat');
zzb =zz;
rrb =interp1(zb,rb,zzb,'spline');
psib =interp1(zb,psib,zzb,'spline');

%% plot
% axes('position',[0.11,0.76,0.85,0.22]);%fig-2 [left bottom width height]
% axes('position',[0.07,0.59,0.425,0.38]);%fig-1[left bottom width height]
axes('position',[0.11,0.75,0.85,0.218]);%fig-1
plot(zb,psib,'m','linewidth',2);
grid on;box on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gca,'xtick',[-2:1:2]);
ylabel('$\psi_w$','interpreter','latex','FontSize',25);
text('Interpreter','latex','String','$(a)$','Position',[-1.4 -1.5],'color','k','FontSize',20)

% axes('position',[0.11,0.53,0.85,0.22]);%fig-2 
% axes('position',[0.555,0.59,0.425,0.38]);%fig-1[left bottom width height]
axes('position',[0.11,0.52,0.85,0.218]);%fig-1
indtmp_1= psi_matlab>0; psio_1=psi_matlab; psio_1(indtmp_1)=NaN;
indtmp_1=find(psi_matlab<0); psii_1=psi_matlab; psii_1(indtmp_1)=NaN;
contour(Y,X,psio_1,15,'linewidth',1); hold on;     %外部区域
contour(Y,X,psii_1,15,'linewidth',1); hold on;   %内部区域  
% contour(Y,X,psi_matlab,[0,0,0],'r--'); hold on;    
set(gca,'xtick',[-2:1:2]);colorbar('east');
ylabel('$r[m]$','interpreter','latex','FontSize',25);
text('Interpreter','latex','String','$(b)$','Position',[-1.4 0.14],'color','k','FontSize',20)

% axes('position',[0.11,0.3,0.85,0.22]);%fig-3 
% axes('position',[0.07,0.11,0.425,0.38]);%fig-1[left bottom width height]
axes('position',[0.11,0.293,0.85,0.218]);%fig-1
indtmp= psi_grass>0; psio=psi_grass; psio(indtmp)=NaN;
indtmp=find(psi_grass<0); psii=psi_grass; psii(indtmp)=NaN;
contour(Y,X,psio,15,'linewidth',1); hold on;     %芯部区域
contour(Y,X,psii,15,'linewidth',1); hold on;    %四周区域
contour(Y,X,psi_grass,[0,0,0],'r--'); hold on;    %四周区域
set(gca,'xtick',[-2:1:2]);colorbar('east');
ylabel('$r[m]$','interpreter','latex','FontSize',25);
text('Interpreter','latex','String','$(c)$','Position',[-1.4 0.14],'color','k','FontSize',20);
xlabel('$z[m]$','interpreter','latex','FontSize',25);

psiin=psi; psiin(psi_matlab<0)=NaN;
% axes('position',[0.11,0.07,0.85,0.22]);%fig-4
% axes('position',[0.555,0.11,0.425,0.38]);%fig-1[left bottom width height]
axes('position',[0.11,0.064,0.85,0.218]);%fig-1
contour(Y,X,psi,13,'linewidth',1); hold on;
% pcolor(Y,X,psi); shading interp; hold on; 
xlabel('$z[m]$','interpreter','latex','FontSize',25);
ylabel('$r[m]$','interpreter','latex','FontSize',25);
h = colorbar('east');
h.TickLabels = (-4:2:4);
h.Label.String = '\times 10^{-5}';
text('Interpreter','latex','String','$(d)$','Position',[-1.4 0.14],'color','k','FontSize',20)
colormap('jet')

