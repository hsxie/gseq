% 19-06-20 14:40

close all;
hf=figure('unit','normalized','Position',[0.4 0.2 0.6 0.5],...
   'DefaultAxesFontSize',14);
% global rs;
% subplot(233);
% surf(Z,R,-Jzeta);
% xlabel('z');ylabel('r');zlabel('-J');

% subplot(234);
% surf(Z,R,P);
% xlabel('z');ylabel('r');zlabel('P');
% Zs=0.566;

Zs=0.5;
subplot(2,3,1:2); 
% Jzetar=-Jzeta(:,floor(Zs*nZ));
ido=find(psi==min(min(psi)));
indtmp= psi>0; psio=psi; psio(indtmp)=NaN;
indtmp=find(psi<0); psii=psi; psii(indtmp)=NaN;
contour(Z,R,psio,7,'linewidth',2); hold on;
contour(Z,R,psii,15,'linewidth',2); hold on;colorbar;
% [c1,h1]=contour(Z,R,psi,[-2E-6,-2E-6],'r--','linewidth',3); hold on; % highligth separatix
[c1,h1]=contour(Z,R,psi,[-0,0],'r--','linewidth',3);
plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3); % O-point
text(0.01*zw,R(ido(1)),'O point','fontsize',13);

ylabel('r (m)');
zx1=c1(1,2:(end-1)); rx1=c1(2,2:(end-1)); zx_no=size(zx1,2);
indx=find(abs(zx1)==min(abs(zx1)));
title(['\psi, I_z=',num2str(Izeta/1e6),'MA, r_o=',num2str(R(ido(1))),...
    ', r_x=',num2str(rx1(indx(1))),', \psi_{min}=',num2str(min(min(psi)),3)]);
% strp='p(\psi)=C*sech[C_2*(\psi-C_3)]';
strp='p(\psi)=-C\cdotexp[-C1\cdot\psi] \times (C2\cdot\psi+1)^2';
xlabel(['z(m), C=',num2str(Cq),...
    ', C1=',num2str(C1q),', C2=',num2str(C2q),...
    ', ',strp]);
ylim([0,rm+0.01]);
% set(gca,'ytick',[0 0.2 0.4 0.6 0.8]);
set(gca,'ytick',[0:0.05:0.20]);

min(psi(:,floor(Zs*nZ)))


% if Zs==0.5
    rs =rx1(indx(1));
    
% else
%     zx1_new=zx1(5:zx_no-5);rx1_new=rx1(5:zx_no-5);
%     rs = interp1(zx1_new,rx1_new,zz(floor(Zs*nZ)),'spline');
% end


% plot([Z(:,floor(Zs*nZ)) Z(:,floor(Zs*nZ))],[0 rrb(floor(Zs*nZ))],'k--','linewidth',3);hold all;

    
plot(zzb,rrb,'k-','linewidth',3);
subplot(233)
R0=floor(R(ido(1))/rw*nR);
yyaxis left;
plot(zz,Bz(R0,:),'-h','LineWidth',1.5,...
'MarkerSize',4);box on;grid on;hold all;
ylabel('$B_{z} \ (T)$','interpreter','latex','FontSize',17);
yyaxis right;
plot(zz,Br(R0,:),'-h','LineWidth',1.5,...
'MarkerSize',4);box on;grid on;hold all;
set(gca,'Fontsize',17);
xlabel('$z\ (m),r=r_{o}$','interpreter','latex','FontSize',17);
ylabel('$B_{r} \ (T)$','interpreter','latex','FontSize',17);
xlim([-zw,zw]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
h=legend('$B_{z}$','$B_{r}$');set(h,'interpreter','latex');
legend('boxoff');

subplot(234);
Jzetar=-Jzeta(:,floor(Zs*nZ));

% tot_I=sum(Jzetar(indx(1)).*dR);
tot_I1=0;
for ii=1:indx(1)
    tot_I1=tot_I1+Jzetar(ii)*dR;
end
% tot_I=tot_I1*L;

L=max(c1(1,:))-min(c1(1,:));
tot_I=tot_I1*L;
Jz1=0.0;
for jz=1:nZ % 19-06-23 22:00 update
    for jr=1:(jrb(jz)-1)
        if (psio(jr,jz)<=0)
            Jz1=Jz1+Jzeta(jr,jz)*dR*dZ;
        end
    end
end
fprintf('I at the midplane inside the separatrix is: %5.4d \n',tot_I);
fprintf('I inside the separatrix is: %5.4d \n',Jz1);


Pr=P(:,floor(Zs*nZ));
J0=max(abs(Jzetar));
P0=max(abs(Pr));
plot(rr,Jzetar/J0,'--',rr,Pr/P0,'linewidth',2);grid on;hold all;
plot([rs rs],[0 1],'-.k','LineWidth',1.5);hold all;
set(gca,'Fontsize',17);
xlabel('$r\ (m), z=0$','interpreter','latex','FontSize',17);
ylabel(['P_0=',num2str(P0,4),', J_0=',num2str(J0,4)]);
h=legend('$-J/J_0$','$P/P_0$');set(h,'interpreter','latex');
legend('boxoff');
xlim([0,0.1]);
set(gca,'xtick',[0:0.02:0.1]);
text('Interpreter','latex','String','$Separatrix$', 'Position',[rs-0.02 0.15],'FontSize',17)
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

P_ave=0;
r_dr=0;
Rend=0.4;
for ii=1:size(rr,2)
    if rr(ii)<=rs
        PP(ii)=Pr(ii);
        rrs(ii)=rr(ii);
        No=ii;
    end
    
%     if rr(ii)<=Rend
%         rrzw(ii)=rr(ii);
%         Bwz(ii)=Bz(ii,floor(0.99*nZ));
%     end
end

for ii=1:size(rr,2)-1
    if rr(ii+1)<=Rend
        rrzw(ii)=rr(ii);
        Bwz(ii)=Bz(ii,floor(0.99*nZ));
    end
end

for jj=No+1:size(rr,2)
    rrw(jj-No)=rr(jj);
    Bw(jj-No)=Bz(jj,floor(0.5*nZ));
end

fprintf('********************************************************* \n');
P_ave=trapz(rrs,rrs.*PP)/trapz(rrs,rrs);
mu0=4*pi*1e-7;
Bx=0.0341;
ave_beta=2*mu0*P_ave/(Bx*Bx);
fprintf('ave_beta = %5.4f \n',ave_beta);

% K=0.530509;
% fprintf('ave_beta = %5.4f \n',tanh(K)/K);

Bz_o=trapz(rrw,rrw.*Bw)/trapz(rrw,rrw);
Bz_w=trapz(rrzw,rrzw.*Bwz)/trapz(rrzw,rrzw);
% rwm=Bz_o*rs*rs/(Bz_o-Bz_w);
% rwm=sqrt(rw);


N_rs=floor(rs/rw*nR);
ps=(Pr(N_rs)+Pr(N_rs+1))/2;
dpdr=(Pr(N_rs+1)-Pr(N_rs))/(rr(N_rs)-rr(N_rs+1));
lp=ps/dpdr;
fprintf('Lp = %5.4f \n',lp);
fprintf('********************************************************* \n');

subplot(235);
yyaxis left;
plot(rr,Bz(:,floor(Zs*nZ)),'-h','LineWidth',1.5,...
'MarkerSize',4);box on;grid on;hold all;
plot([rs rs],[1.2*min(Bz(:,floor(Zs*nZ))) max(Bz(:,floor(Zs*nZ)))]...
    ,'-.k','LineWidth',1.5);hold all;
ylabel('$B_{z} \ (T)$','interpreter','latex','FontSize',17);
yyaxis right;
plot(rr,Br(:,floor(Zs*nZ)),'-h','LineWidth',1.5,...
'MarkerSize',4);box on;grid on;hold all;
set(gca,'Fontsize',17);
xlabel('$r\ (m),z=0$','interpreter','latex','FontSize',17);
ylabel('$B_{r} \ (T)$','interpreter','latex','FontSize',17);
xlim([0,0.1]);
set(gca,'xtick',[0:0.02:0.1]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
% legend('B_z','B_r');legend('boxoff');
h=legend('$B_{z}$','$B_{r}$');set(h,'interpreter','latex');
legend('boxoff');

% text('Interpreter','latex','String','$Separatrix$', 'Position',...
%     [rs-0.1 0.1*max(Br(:,floor(Zs*nZ)))],'FontSize',17)

subplot(236);
yyaxis left;
Z_end=0.99;
plot(rr,Bz(:,floor(Z_end*nZ)),'-h','LineWidth',1.5,...
'MarkerSize',4);box on;grid on;hold all;
ylabel('$B_{z} \ (T)$','interpreter','latex','FontSize',17);
text('Interpreter','latex','String','$Separatrix$', 'Position',...
    [rs-0.02 1.1*min(Bz(:,floor(Z_end*nZ)))],'FontSize',17)
yyaxis right;
plot(rr,Br(:,floor(Z_end*nZ)),'-h','LineWidth',1.5,...
'MarkerSize',4);box on;grid on;hold all;
set(gca,'Fontsize',17);
xlabel('$r\ (m),z=0.99 z_{w}$','interpreter','latex','FontSize',17);
ylabel('$B_{r} \ (T)$','interpreter','latex','FontSize',17);
xlim([0,0.15]);
set(gca,'xtick',[0:0.05:0.15]);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
h=legend('$B_{z}$','$B_{r}$');set(h,'interpreter','latex');
legend('boxoff');



set(gcf, 'PaperPositionMode', 'auto');
% RemovePlotWhiteArea(gca);

