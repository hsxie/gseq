%contour or pcolor data from grass_ft
close all;clear;clc;

figure('position',[50 100 600 300]);
set(gca,'nextplot','replacechildren');
set(gcf,'DefaultAxesFontSize',15);

rmin =0;
zmin =-1.5;
rmax =1.0;
zmax =1.5;

data =importdata('data/000000_pse.DAT');
r =linspace(rmin,rmax,size(data,1));
z =linspace(zmin,zmax,size(data,2));
[X,Y]=ndgrid(r,z);

rr =linspace(rmin,rmax,size(data,1));
zz =linspace(zmin,zmax,size(data,2));
[X1,Y1] =ndgrid(rr,zz);

for i=1:50:2001
    name2=num2str(i-1,'%06d');
    xm =['data/',name2,'_pse.DAT'];
    data =importdata(xm);
    
    Z1 =griddata(X,Y,data,X1,Y1);
    set(gca,'Fontsize',15);
    contour(Y1,X1,Z1,30,'linewidth',1);     %Ð¾²¿ÇøÓò
    title(['FRC G-S equilibrium ', ' it=',num2str(i-1)]);
    box on;grid on;colorbar;caxis([-4 0]);
    xlim ([zmin,zmax]);ylim ([0,1]);
    set(gca,'XTick',[zmin:0.5:zmax]);
    set(gca,'YTick',[0:0.2:1]);
    xlabel('$Z$','interpreter','latex','FontSize',17);
    ylabel('$R$','interpreter','latex','FontSize',17);
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);    
    pause(0.1);
    
    frame=getframe(gcf);
    imind=frame2im(frame);
    [imind,cm] = rgb2ind(imind,256);
    if i==1
         imwrite(imind,cm,'test.gif','gif', 'Loopcount',inf,'DelayTime',1e-4);
    else
         imwrite(imind,cm,'test.gif','gif','WriteMode','append','DelayTime',1e-4);
    end
end
% imwrite(imind,cm,'test.gif','GIF', 'Loopcount',inf,'DelayTime',1);

