%contour or pcolor data from grass_ft
close all;clear;clc;

figure('position',[50 100 600 300]);
set(gca,'nextplot','replacechildren');
set(gcf,'DefaultAxesFontSize',15);

rmin =0;
zmin =-1.5;
rmax =1.0;
zmax =1.5;

nplot =1;
if nplot ==1
    data =importdata('data/ps_t.dat');
elseif nplot ==2
    data =importdata('data/ps_pl.dat');
elseif nplot ==3
    data =importdata('data/ps_t_ref.dat');
end 

r =linspace(rmin,rmax,size(data,1));
z =linspace(zmin,zmax,size(data,2));
[X,Y]=ndgrid(r,z);

rr =linspace(rmin,rmax,size(data,1));
zz =linspace(zmin,zmax,size(data,2));
[X1,Y1] =ndgrid(rr,zz);

set(gca,'Fontsize',15);
Z1 =griddata(X,Y,data,X1,Y1);
% contour(Y1,X1,Z1,200);

indtmp= Z1>0; psio=Z1; psio(indtmp)=NaN; 
indtmp=find(Z1<0); psii=Z1; psii(indtmp)=NaN;
contour(Y1,X1,psio,15,'linewidth',1); hold on;     %芯部区域
contour(Y1,X1,psii,15,'linewidth',1); hold on;    %四周区域
contour(Y1,X1,Z1,[0,0,0],'r--'); hold on;    %四周区域

hold all;box on;grid on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
shading flat;colorbar;
xlim ([zmin,zmax]);ylim ([0,1.01]);
set(gca,'XTick',[zmin:0.5:zmax]);
xlabel('$Z$','interpreter','latex','FontSize',17);
ylabel('$R$','interpreter','latex','FontSize',17);

if nplot ==1
    title('$Grad-Shafranov \ eq \ \psi_{e}$','interpreter','latex','FontSize',17);
elseif nplot ==2
    title('$Grad-Shafranov \ eq \ \psi_{p}$','interpreter','latex','FontSize',17);
elseif nplot ==3
    title('$Grad-Shafranov \ eq \ \psi_{tot}$','interpreter','latex','FontSize',17);
end 

psiq=data;
filename='psi_grass';
save([filename,'.mat'],'Z1');

% set(gca,'Fontsize',15);
% plot(z,data(129,:),'m');hold all;box on;grid on;
% plot(linspace(-0.25,-0.25,100),linspace(-4,-1,100),'r--');hold all;
% plot(linspace(0.25,0.25,100),linspace(-4,-1,100),'r--');hold all;
% 
% plot(linspace(-0.5,-0.5,100),linspace(-4,-1,100),'b--');hold all;
% plot(linspace(0.5,0.5,100),linspace(-4,-1,100),'b--');hold all;
% 
% xlim ([zmin,zmax]);ylim ([-4,-1]);
% set(gca,'XTick',[zmin:0.5:zmax]);
% set(gca,'YTick',[-4:1:-1]);
% xlabel('$Z$','interpreter','latex','FontSize',17);
% ylabel('$\psi_{e}$','interpreter','latex','FontSize',17);
% set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
% 
% text('Interpreter','latex','String','$r =r_{wall}$', ...
% 'Position',[0,-2.5],'FontSize',15) ;
% 
% 
max(max(data))
