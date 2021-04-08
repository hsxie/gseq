% 2020-07-31 15:08, Hua-sheng Xie, ENN
% Solve 3PE parameters, (betas,ls)->(alpha,sigma,psim)
% 20-08-02 08:02 test seems ok
close all;clear;clc;

global betas ls alpha pm tol alpha0 ls0 psim0 psig;
pm=1;
tol=1e-8;

betas=0.6;
% betas=0.6655;
alpha0=8*log((sqrt(1-betas)+1)/sqrt(betas));
ls0=1/(alpha0*sqrt(1-betas));
psim0=log(betas)/alpha0;
% ls=0.75;
% ls=0.00823/0.0561;
ls=0.1;
% ls=ls0;

% alpha=2.082049233064903;
% sigma=1.687721726802611;
% psim=-0.056010101372585;
% ls=0.147634407671656;
% betas=0.6655;

% alpha=1/(ls*sqrt(1-betas));

psig=-0.035;

options = optimset('TolFun',tol,'Tolx',tol,'Maxfuneval',10000,'MaxIter',5000);

x0=0.8;
x=fsolve(@f3pe,x0,options)
F=f3pe(x)
sigma=real(x)

if(1==1)
    aa=-5.0:0.11:3.1;ff=0.*aa;
    for ja=1:length(aa)
        ff(ja)=f3pe(aa(ja));
    end
    plot(aa,real(ff),aa,imag(ff),'--','linewidth',2);
    xlabel('\sigma');ylabel('F_3(\sigma)');
    legend('Re[F]','Im[F]','location','best');legend('boxoff');
    title(['\beta_s=',num2str(betas),', l_s=',num2str(ls),...
        ', l_{s0}^{RR}=',num2str(ls0),', pm=',num2str(pm)]);

end

if(1==0)
    sigma=1;
alpha=((16*sigma*ls*sqrt(1-betas)+7)+pm*sqrt((16*sigma*ls*sqrt(1-betas)+7)^2-...
    4*sigma*48*ls*sqrt(1-betas)))/(2*48*ls*sqrt(1-betas));

aa=-0.3:0.01:0;
pp=(betas-8*betas*alpha^2/(7*alpha-sigma)*(8*aa.*(1+4*sigma*aa)+...
    (sigma+alpha)*(1-exp(16*alpha*aa))/(8*alpha^2)));


plot(aa,real(pp),aa,imag(pp));

end

alpha=((16*sigma*ls*sqrt(1-betas)+7)+pm*sqrt((16*sigma*ls*sqrt(1-betas)+7)^2-...
    4*sigma*48*ls*sqrt(1-betas)))/(2*48*ls*sqrt(1-betas));

F1=@(psim)(betas-8*betas*alpha^2/(7*alpha-sigma)*(8*psim.*(1+4*sigma*psim)+...
    (sigma+alpha)*(1-exp(16*alpha*psim))/(8*alpha^2)))-1+1e-6;
% F11=@(psim)(F1(psim)).^2;
% pp=0.0:-0.01:-0.3;plot(pp,F1(pp)+1/betas);
%

% options2 = optimset('TolFun',1e-8,'Tolx',1e-8,'Display','iter');
options2 = optimset('TolFun',tol,'Tolx',tol,'Display','off');
psim=fsolve(F1,psig,options2);
% psim=fzero(F1,psig,options2);
% psim=fzero(F1,-0.005,options2);
% psim=fzero(F1,psim0,options2);
% psim=fminsearch(F11,psig,options2);
F1(psim)

% F=f3pe(x)
% pp=0:-0.01:-0.3;plot(pp,F1(pp));

%%
close all;

fP=@(psi)betas-8*betas*alpha^2/(7*alpha-sigma).*(...
    (psi<0).*(8*psi.*(1+4*sigma.*psi)+...
    (sigma+alpha).*((1-exp(16*alpha.*psi))/(8*alpha^2)))+...
    (psi>=0).*((1-exp(-8*alpha.*psi))/(alpha)-...
    (sigma+alpha).*(1-exp(-16*alpha.*psi))/(8*alpha^2) ));

fdP=@(psi)-8*betas*alpha^2/(7*alpha-sigma)*(...
    (psi<0).*(8*(1+8*sigma*psi)+...
    (sigma+alpha)*(-2*exp(16*alpha*psi))/(alpha))+...
    (psi>=0).*(8*exp(-8*alpha*psi)+...
    (sigma+alpha)*(-2*exp(-16*alpha*psi))/(alpha)));

fpsi=@(psi)1./sqrt(1.0-fP(psi)+0e-4);

dr=0.001; rmax=1.5;
rr=0:dr:rmax;
% nr=rmax/dr;
nr=length(rr);
psi=0.*rr;

for jr=2:nr
    r=rr(jr);
    if(r<=1/sqrt(2))
        FF1=@(psi)integral(fpsi,0,psi,'RelTol',tol,'AbsTol',tol)+r^2/2;
%         psi(jr)=fsolve(FF1,-0.01,options2);
        psi(jr)=fsolve(FF1,real(psi(jr-1)),options2);
    else
        FF2=@(psi)integral(fpsi,psim,psi,'RelTol',tol,'AbsTol',tol)-(r^2-1/2)/2;
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

% to update r>rs
Pr=fP(psi);

dP=fdP(psi);
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
title(['l_{s0}=',num2str(ls0),', \alpha=',num2str(alpha),', pm=',num2str(pm)]);
subplot(223);
plot(rr,Jt,'linewidth',2);
xlabel('r');ylabel('J_\theta');
title(['\sigma=',num2str(sigma),', h=',num2str(h),', \psi_g=',num2str(psig)]);
subplot(224);
plot(rr,Bz,'linewidth',2);
xlabel('r');ylabel('B_z');
title(['\psi_m=',num2str(psim),', \alpha_0=',num2str(alpha0),', x_0=',num2str(x0)]);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);
print(gcf,'-dpng',['3PE_betas=',num2str(betas),'_ls=',num2str(ls),...
    '_pm=',num2str(pm),'_alpha=',num2str(alpha),...
    '_tol=',num2str(tol),'.png']);

%%
if(1==0)
alpha=2.082049233064903;
sigma=1.687721726802611;
psim=-0.056010101372585;
% ls=0.008230983716240/0.0561;
% alpha=0.99;
% sigma=0.8023;
% psim=-0.117793429722030;
ls=0.147634407671656;
betas=0.6655;

psi=psim:0.001:4*abs(psim);
plot(psi,fP(psi));

F1=@(psim)fP(psim)-1;

% alpha1=((16*sigma*ls*sqrt(1-betas)-7)+pm*sqrt((16*sigma*ls*sqrt(1-betas)-7)^2+...
%     4*sigma*48*ls*sqrt(1-betas)))/(2*48*ls*sqrt(1-betas))
alpha1=((16*sigma*ls*sqrt(1-betas)+7)+pm*sqrt((16*sigma*ls*sqrt(1-betas)+7)^2-...
    4*sigma*48*ls*sqrt(1-betas)))/(2*48*ls*sqrt(1-betas))

Falpha=ls*16*alpha*(3*alpha-sigma)*sqrt(1-betas)-(7*alpha-sigma)
Falpha2=1/(8*ls)-8*alpha^2/(7*alpha-sigma)*(1-(sigma+alpha)/(4*alpha))*sqrt(1-betas)
F1(psim)
end

datq_3pe=[rr;Pr;psi;Bz;Jt].';
save('3pe_datq.mat','datq_3pe','betas','ls');
