% 2020-07-28 21:05, Hua-sheng Xie, ENN
% Solve MRR parameters, (betas,ls)->(alpha,sigma,psim)
% 20-07-29 16:51 update
% 20-07-31 10:17 modified to MRR1
close all;clear;clc;

global betas ls alpha n tol alpha0 ls0 psim0 psig;
n=2;
tol=1e-8;

betas=0.6655;
% betas=0.6;
alpha0=8*log((sqrt(1-betas)+1)/sqrt(betas));
ls0=1/(alpha0*sqrt(1-betas));
psim0=log(betas)/alpha0;
ls=0.148;
% ls=ls0;

alpha=1/(ls*sqrt(1-betas));

psig=-0.05;
% psig=-0.06;

options = optimset('TolFun',tol,'Tolx',tol,'Maxfuneval',10000,'MaxIter',5000);

x0=71;
% x0=3.0;
x=fsolve(@fmrr2,x0,options)
F=fmrr2(x)
sigma=real(x)
F1=@(psim)exp(-alpha*psim-sigma*psim.^n)-1/betas;

% options2 = optimset('TolFun',1e-8,'Tolx',1e-8,'Display','iter');
options2 = optimset('TolFun',tol,'Tolx',tol,'Display','off');
psim=fsolve(F1,psig,options2);
% psim=fzero(F1,-0.005,options2);
% psim=fzero(F1,psim0,options2);
F1(psim)

if(1==1)
    aa=-50.0:1.5:100.1;ff=0.*aa;
    for ja=1:length(aa)
        ff(ja)=fmrr2(aa(ja));
    end
    plot(aa,real(ff),aa,imag(ff),'--','linewidth',2);
    xlabel('\sigma');ylabel('F_3(\sigma)');
    legend('Re[F]','Im[F]','location','best');legend('boxoff');
    title(['\beta_s=',num2str(betas),', l_s=',num2str(ls),...
        ', l_{s0}^{RR}=',num2str(ls0)]);

end

% aa=-10:0.1:50;ff=0.*aa;
% for ja=1:length(aa)
%     ff(ja)=fmrr(aa(ja));
%     
% end
% plot(aa,ff);

% F=fmrr(x)
% pp=0:-0.01:-0.3;plot(pp,F1(pp));

%%
close all;

fpsi=@(psi)1./sqrt(1.0-betas*exp(-alpha*psi-sigma*psi.^n));

dr=0.002; rmax=1.5/1;
rr=0:dr:rmax;
% nr=rmax/dr;
nr=length(rr);
psi=0.*rr;

for jr=2:nr
    r=rr(jr);
    if(r<=1/sqrt(2))
        FF1=@(psi)integral(fpsi,0,psi,'RelTol',1e-8,'AbsTol',1e-8)+r^2/2;
%         psi(jr)=fsolve(FF1,-0.01,options2);
        psi(jr)=fsolve(FF1,real(psi(jr-1)),options2);
    else
        FF2=@(psi)integral(fpsi,psim,psi,'RelTol',1e-8,'AbsTol',1e-8)-(r^2-1/2)/2;
%         psi(jr)=fsolve(FF2,psim0,options2);
        psi(jr)=fsolve(FF2,0.99*real(psi(jr-1)),options2);
    end
end

%%
% psi(abs(imag(psi))>=1e-6)=NaN;
close all;

figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.02,0.2,0.5,0.6]);
fP=@(psi)betas*exp(-alpha*psi-sigma*psi.^n);
Pr=betas*exp(-alpha*psi-sigma*psi.^n);
fdP=@(psi)-betas*exp(-alpha*psi-sigma*psi.^n).*(alpha+n*sigma*psi.^(n-1));
dP=-betas*exp(-alpha*psi-sigma*psi.^n).*(alpha+n*sigma*psi.^(n-1));
Jt=-rr.*dP;
Bz=-sqrt(1-Pr).*(rr<=1/sqrt(2))+sqrt(1-Pr).*(rr>1/sqrt(2));

ids=find(abs(rr-1)==min(abs(rr-1)));ids=ids(1);
lsp=real(-(Pr(ids)/((Pr(ids+1)-Pr(ids-1))/(rr(ids+1)-rr(ids-1)))));

ind=find(rr<=1.0);
h=real(fdP(psim)/sum(2*dP(ind).*rr(ind)*dr));

subplot(221);
plot(rr,psi,'linewidth',2);
xlabel('r');ylabel('\psi');
title(['\beta_s=',num2str(betas),', l_s=',num2str(ls),...
    ', l_{sp}=',num2str(lsp,3)]);
subplot(222);
plot(rr,Pr,'linewidth',2);
xlabel('r');ylabel('P');
title(['l_{s0}=',num2str(ls0),', \alpha=',num2str(alpha),', \psi_g=',num2str(psig)]);
subplot(223);
plot(rr,Jt,'linewidth',2);
xlabel('r');ylabel('J_\theta');
title(['\sigma=',num2str(sigma),', h=',num2str(h),', n=',num2str(n)]);
subplot(224);
plot(rr,Bz,'linewidth',2);
xlabel('r');ylabel('B_z');
title(['\psi_m=',num2str(psim),', \alpha_0=',num2str(alpha0),', x_0=',num2str(x0)]);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);
print(gcf,'-dpng',['MRR2_betas=',num2str(betas),'_ls=',num2str(ls),...
    '_n=',num2str(n),'_alpha=',num2str(alpha),...
    '_tol=',num2str(tol),'.png']);

datq_mrr2=[rr;Pr;psi;Bz;Jt].';
% save('mrr2_datq.mat','datq_mrr2','betas','ls');
% save('mrr2_datq.mat','datq_mrr2','betas','ls');