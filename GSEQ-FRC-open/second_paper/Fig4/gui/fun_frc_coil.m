% Hua-sheng Xie, huashengxie@gmail.com, CCF-ENN, 2019-06-13 11:00
% Draw FRC coil
function [runtime,psi0,Br0,Bz0]=fun_frc_coil(r0,z0,rw,zw,rl,zl,hl,wl,Il,Nl,imethod,hAxes1,hAxes2,hAxes3)

disp('Running fun_frc_coil.m ...');
runtime=cputime;
% imethod=3; % =3, YuanBS book p110; =1, Green function; =2, nr*nz Green function
% r0
% z0
% rl
% zl
% hl
% wl
% Il

[psi0,Br0,Bz0]=fun_Bfield(r0,z0,rl,zl,hl,wl,Il,Nl,imethod);



set(gcf,'currentaxes',hAxes1);
cla;
%   rw=1.0; zw=2.5;

%   rl=[0.7,1.1,0.7]*rw; % r position of coils
%   zl=[-1.15,0.0,1.15]*zw; % z position of coils
%   hl=[0.15,0.1,0.15]*rw; % r position of coils
%   wl=[0.1,1.6,0.1]*zw; % z position of coils
%   Il=[10e5,1.0e5,10e5]; % current of coils, Ampere
  

nl=length(rl); % # of coils
if((nl~=length(zl)) ||(nl~=length(wl)) ||(nl~=length(hl)) ||(nl~=length(Il)) ||(nl~=length(Nl)) )
    warndlg('Length of rl, zl, hl, wl, Il, Nl not consistient!');
    return;
end

% fontsize=25;
% figure('unit','normalized','DefaultAxesFontSize',fontsize,...
%     'DefaultAxesFontWeight','bold',...
%     'DefaultAxesLineWidth',2,'defaulttextinterpreter','latex',...
%     'position',[0.1,0.2,0.28,0.25],'color', [1, 1, 1]);

nR=2^5+1; nZ=2^6+1; % grid numbers
% nR=2^3+1; nZ=2^4+1; % grid numbers

Rlow=0; Rup=rw; Zlow=-zw; Zup=zw;
rr=linspace(Rlow,Rup,nR); zz=linspace(Zlow,Zup,nZ);
% dR=rr(2)-rr(1); dZ=zz(2)-zz(1);
[R,Z]=ndgrid(rr,zz); %dR2=dR*dR; dZ2=dZ*dZ;
psi=0.*R;
for jr=1:nR
    for jz=1:nZ
        [psi(jr,jz),~,~]=fun_Bfield(rr(jr),zz(jz),rl,zl,hl,wl,Il,Nl,imethod);
    end
end
psi(abs(psi)>1e5)=NaN;
% contour(Z,R,psiw_coil);
% pcolor(Z,R,real(psi)); shading interp; hold on; 
contour(Z,R,real(psi),100);hold on;

contour(Z,-R,real(psi),100); hold on;colorbar('east');


for jl=1:nl
rectangle('position',[zl(jl)-0.5*wl(jl),rl(jl)-0.5*hl(jl),...
    wl(jl),hl(jl)],'FaceColor', [0.5 0.5 0.5],'linewidth',1);
hold on;
rectangle('position',[zl(jl)-0.5*wl(jl),-rl(jl)-0.5*hl(jl),...
    wl(jl),hl(jl)],'FaceColor', [0.5 0.5 0.5],'linewidth',1);
hold on;
% if(Il(jl)>0)
%     plot(zl(jl),rl(jl),'r.',zl(jl),-rl(jl),'rx','linewidth',2);
% else
%     plot(zl(jl),rl(jl),'rx',zl(jl),-rl(jl),'r.','linewidth',2);
% end
zw=0.75;rw=0.15;
zz=-zw:0.01*zw:zw;
rr=0.14+0.*zz;
rx=0.*zz;
plot(zz,rr,'r-','linewidth',2); hold on;
plot(zz,-rr,'r-','linewidth',2);hold on;
plot(zz,rx,'k--','linewidth',2);hold on;
plot([-0.4301,0.4301],[0,0],'r--','linewidth',4);hold on;
set(gca,'ytick',[-0.15:0.1:0.15]);set(gca,'xtick',[-0.75:0.25:0.75]);
end

xlim([Zlow,Zup]);
ylim([-Rup,Rup]);
xlabel('$z (m)$','interpreter','latex','FontSize',25);
ylabel('$r (m)$','interpreter','latex','FontSize',25);
% axis equal; 
box on;


set(gcf,'currentaxes',hAxes2);
cla;

Br1=0.*rr;
Bz1=0.*rr;
for jr=1:nR
    [~,Br1(jr),Bz1(jr)]=fun_Bfield(rr(jr),z0,rl,zl,hl,wl,Il,Nl,imethod);
end
B1=sqrt(Br1.^2+Bz1.^2);

plot(rr,Br1,'.',rr,Bz1,'--',rr,B1,':','linewidth',2);
xlabel(['r (m), z=',num2str(z0),'m']);ylabel('B (T)');
legend('B_r','B_z','B');legend('boxoff');

% ylim([-5e-3,1.5e-2]);
xlim([0,rw]);

set(gcf,'currentaxes',hAxes3);
cla;

Br2=0.*zz;
Bz2=0.*zz;
for jz=1:nZ
    [~,Br2(jz),Bz2(jz)]=fun_Bfield(r0,zz(jz),rl,zl,hl,wl,Il,Nl,imethod);
end
B2=sqrt(Br2.^2+Bz2.^2);

plot(zz,Br2,'.',zz,Bz2,'--',zz,B2,':','linewidth',2);
legend('B_r','B_z','B');legend('boxoff');
xlabel(['z (m), r=',num2str(r0),'m']); ylabel('B (T)');
% ylim([-0.05,0.25]);
xlim([Zlow,Zup]);

runtime=cputime-runtime;

set(gcf,'currentaxes',hAxes1);
title(['method=',num2str(imethod),', runtime=',num2str(runtime),'s']);

disp('Finish run fun_frc_coil.m.');

end