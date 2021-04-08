% 2020-08-01 12:27, Hua-sheng Xie, ENN
% Solve 2PE parameters, (betas,ls)->(sigma,psim)
% 20-08-04 00:04 update to Steinhauer2009, with r>rs
close all;clear;clc;

global betas ls alpha0 ls0 psim0;

betas=0.6655;
alpha0=8*log((sqrt(1-betas)+1)/sqrt(betas));
ls0=1/(alpha0*sqrt(1-betas));
psim0=log(betas)/alpha0;
ls=0.148;

sigma=(betas/(4*(1-betas)*ls)-1)/4;
bs=sqrt(1-betas);
sigma2=(1+4*sigma)*bs/(1-bs);
% psim=sqrt(1-betas)/(8*sigma)*(exp(-sigma)-1);

fb=@(u)(abs(u)<=1).*(bs.*u.*exp(sigma*(u.^4-1)))+...
    (u>1).*(1-(1-bs)*exp(-sigma2*(u-1)));
fp=@(u)1-(fb(u)).^2;
dbdu=@(u)(abs(u)<=1).*(bs.*exp(sigma*(u.^4-1)).*(1+4*sigma*u.^4))+...
    (u>1).*((1-bs)*sigma2*exp(-sigma2*(u-1)));

dr=0.005;
rr=0:dr:1.5;
u=2*rr.^2-1;
Bz=fb(u);
psi=cumsum(rr.*Bz*dr);
Pr=fp(u);
Jt=4*rr.*dbdu(u); % to check

psim=min(psi);

ids=find(abs(rr-1)==min(abs(rr-1)));ids=ids(1);
lsp=real(-(Pr(ids)/((Pr(ids+1)-Pr(ids-1))/(rr(ids+1)-rr(ids-1)))));

ind=find(rr<=1.0);
h=real(dbdu(0)/sum(2*dbdu(u(ind))*dr));

close all;

figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.02,0.2,0.5,0.6]);

subplot(221);
plot(rr,psi,'linewidth',2);
xlabel('r');ylabel('\psi');
title(['\beta_s=',num2str(betas),', l_s=',num2str(ls),...
    ', l_{sp}=',num2str(lsp,3)]);
subplot(222);
plot(rr,Pr,'linewidth',2);ylim([0,1]);
xlabel('r');ylabel('P');
title(['l_{s0}=',num2str(ls0)]);
subplot(223);
plot(rr,Jt,'linewidth',2);
xlabel('r');ylabel('J_\theta');
title(['\sigma=',num2str(sigma),', h=',num2str(h)]);
subplot(224);
plot(rr,Bz,'linewidth',2);
xlabel('r');ylabel('B_z');
title(['\psi_m=',num2str(psim),', \alpha_0=',num2str(alpha0)]);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);
print(gcf,'-dpng',['2PE09_betas=',num2str(betas),'_ls=',num2str(ls),...
    '_sigma=',num2str(sigma),'_psim=',num2str(psim),'.png']);

datq_2pe=[rr;Pr;psi;Bz;Jt].';
save('2pe_datq_i09.mat','datq_2pe','betas','ls');

