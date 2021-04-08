% 19-06-20 14:40

close all;
fontsize=25;
figure('unit','normalized','DefaultAxesFontSize',fontsize,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,'defaulttextinterpreter','latex',...
    'position',[0.45,0.4,0.5,0.55],'color', [1, 1, 1]);
Bz_plasma=Bz;
psi_plasma=psi;
%%
load('MRR1d.mat');
Be0=Be;
rs0=Rs;
mu0=4*pi*1e-7;
%%
% axes('position',[0.11,0.75,0.85,0.23]);%fig-1
axes('position',[0.065,0.69,0.91,0.3]);%fig-1[left bottom width height]
ido=find(psi==min(min(psi)));
indtmp= psi>0; psio=psi; psio(indtmp)=NaN;
indtmp=find(psi<0); psii=psi; psii(indtmp)=NaN;
contour(Z,R,psio,12,'linewidth',2); hold on;
contour(Z,R,psii,15,'linewidth',2); hold on;colorbar('east');%colormap(hot);
contour(Z,-R,psii,15,'linewidth',2); hold on;
contour(Z,-R,psio,12,'linewidth',2); hold on;
[c1,h1]=contour(Z,R,psi,[-0,0],'r--','linewidth',3);hold on;
contour(Z,-R,psi,[-0,0],'r--','linewidth',3);
plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3);hold on; % O-point
plot(Z(ido(1)),-R(ido(1)),'mx','linewidth',3); % O-point
plot(-Z(ido(1)),R(ido(1)),'mx','linewidth',3);hold on; % O-point
plot(-Z(ido(1)),-R(ido(1)),'mx','linewidth',3); % O-point
set(gca,'Fontsize',fontsize);

xlabel('$z [m]$','interpreter','latex','FontSize',25);
ylabel('$r [m]$','interpreter','latex','FontSize',fontsize);
zx1=c1(1,2:(end-1)); rx1=c1(2,2:(end-1)); zx_no=size(zx1,2);
indx=find(abs(zx1)==min(abs(zx1)));
min(psi(:,Z0))
rs=rs0;

for jl=1:size(zl,2)
    rectangle('position',[zl(jl)-0.5*wl(jl),rl(jl)-0.5*hl(jl),wl(jl),hl(jl)],...
        'FaceColor', [0.5 0.5 0.5],'linewidth',1);hold on;
    rectangle('position',[zl(jl)-0.5*wl(jl),-rl(jl)-0.5*hl(jl),wl(jl),hl(jl)],...
        'FaceColor', [0.5 0.5 0.5],'linewidth',1);hold on;    
end

%plot_chamber
rwc=1;zzc=-zw:0.01*zw:zw;
a=0.15; zmc=0.6*zw; rac=0.3*rwc; rbc=0.4*rwc;
ind1=find(zzc<0);
ind2=find(zzc>=0);
rrc(ind1)=0.5*(rbc-rac)*(tanh(5*(zzc(ind1)+zmc))+1.0)+rac;
rrc(ind2)=0.5*(rbc-rac)*(tanh(5*(zmc-zzc(ind2)))+1.0)+rac;
plot(zzc,rrc,'k-','linewidth',2); hold on;
plot(zzc,-rrc,'k-','linewidth',2);
text('Interpreter','latex','String','$(a)$','Position',[-3.4 0.5],'color','k','FontSize',fontsize)

axes('position',[0.065,0.065,0.425,0.24]);%fig-4[left bottom width height]
Jzetar=-Jzeta(:,jN0);
Pr=P(:,jN0);
J0=max(abs(Jzetar));
P0=max(abs(Pr));
plot(rr*100,Pr,'m-',rr1d.*Rs*100,pr1d.*(Be^2/(2*mu0)),'--','linewidth',3);grid on;hold all;
plot([100*rs 100*rs],1.1*[0 max(Pr)],'-.k','LineWidth',1.5);hold all;
set(gca,'Fontsize',fontsize);
xlabel('$r\ [cm]$','interpreter','latex','FontSize',fontsize);
ylabel('$P$','interpreter','latex','FontSize',fontsize);
h=legend('$2D$','$1D$');set(h,'interpreter','latex');legend('boxoff');
xlim([0,45]);set(gca,'xtick',[0:10:45]);set(gca,'ytick',[0:2.5e4:1e5]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$(d)$','Position',[1.5 9.5e4],'color','k','FontSize',fontsize);
text('Interpreter','latex','String','$Separatrix$', 'Position',...
    [100*rs-8 2e4],'FontSize',20);


% axes('position',[0.59,0.07,0.37,0.26]); %fig-5
axes('position',[0.55,0.065,0.425,0.24]);%fig-5[left bottom width height]
Jzetar=-Jzeta(:,jN0);
plot(rr.*100,-Jzetar,'m-',rr1d.*Rs*100,j1d.*(Be/(2*mu0*Rs)),'--','linewidth',3);grid on;hold all;
plot([100*rs 100*rs],1.3*[0 J0],'-.k','LineWidth',1.5);hold all;
set(gca,'Fontsize',fontsize);
xlabel('$r\ [cm]$','interpreter','latex','FontSize',fontsize);
ylabel('$J_{\theta}$','interpreter','latex','FontSize',fontsize);
h=legend('$2D$','$1D$');set(h,'interpreter','latex');legend('boxoff');
xlim([0,45]);set(gca,'xtick',[0:10:45]);ylim([0,9e6]);set(gca,'ytick',[0:3e6:9e6]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
text('Interpreter','latex','String','$(e)$','Position',[1.5 8e6],'color','k','FontSize',fontsize)


load('psib_FIX.mat'); %read data
rrb=interp1(zb,rb,Z(1,:),'spline');
w_No=fix(rrb./(R(2,1)-R(1,1)));
for i=1:size(w_No,2)
    psi_wn(i)=psi_plasma(w_No(i),i);
    Bz_wn(i)=Bz_plasma(w_No(i),i);
end
% axes('position',[0.12,0.41,0.37,0.26]); %fig-2 [left bottom width height]
axes('position',[0.065,0.39,0.425,0.24]);%fig-2[left bottom width height]
% plot(zb,psib,'r-',Z(1,:),psi_plasma,'g--','linewidth',3);grid on;hold all;
plot(zb,psib,'r-',Z(1,:),psi_wn,'g--','linewidth',3);grid on;hold all;
set(gca,'Fontsize',fontsize);
xlabel('$z [m]$','interpreter','latex','FontSize',fontsize);
ylabel('$\psi \ at \ wall$','interpreter','latex','FontSize',fontsize);
xlim([-3.55,3.5]);set(gca,'xtick',[-3:1:3]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
set(gcf, 'PaperPositionMode', 'auto');
h=legend('$\mathbf{w/o \ plasma}$','$\mathbf{w \ plasma}$');set(h,'interpreter','latex');
legend('boxoff');
text('Interpreter','latex','String','$(b)$','Position',[-3.3 0.038],'color','k','FontSize',fontsize)

% axes('position',[0.59,0.41,0.37,0.26]); %fig-3
axes('position',[0.55,0.39,0.425,0.24]);%fig-3[left bottom width height]
plot(zb,Bz,'r-',Z(1,:),Bz_wn,'g--','linewidth',3);grid on;hold all;
% plot(zb,Bz,'r-',Z(1,:),Bz_plasma,'g--','linewidth',3);grid on;hold all;
ylabel('$B_{z} \ at \ wall$','interpreter','latex','FontSize',fontsize);
xlabel('$z [m]$','interpreter','latex','FontSize',fontsize);
xlim([-3.55,3.5]);set(gca,'xtick',[-3:1:3]);ylim([0,1.5]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.25);
h=legend('$\mathbf{w/o \ plasma}$','$\mathbf{w \ plasma}$');set(h,'interpreter','latex');
legend('boxoff');
text('Interpreter','latex','String','$(c)$','Position',[-3.3 1.35],'color','k','FontSize',fontsize);




