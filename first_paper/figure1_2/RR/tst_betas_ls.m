% 20-07-30 04:52
close all;clear;clc;

global betas ls n tol alpha0 ls0 psim0;
n=2;
tol=1e-8;

bb=0.02:0.02:0.95;
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

figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.02,0.6,0.8,0.3]);

% plot(bb,log10(ll),'linewidth',2);
subplot(131);
semilogy(bb,ll,'linewidth',2);
grid on; box on;
xlabel('\beta_s');ylabel('l_s');

subplot(132);
plot(bb,aa,'linewidth',2);
grid on; box on;
xlabel('\beta_s');ylabel('\alpha');

subplot(133);
plot(bb,pp,'linewidth',2);
grid on; box on;
xlabel('\beta_s');ylabel('\psi_m');
