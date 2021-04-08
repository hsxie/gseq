% 2020-07-28 21:05, Hua-sheng Xie, ENN
% Solve MRR parameters, (betas,ls)->(alpha,sigma,psim)
% 20-07-29 16:51 update
% 20-07-31 10:17 modified to MRR1
% 20-08-02 21:49 modified to RR
close all;clear;clc;

global betas ls alpha;
tol=1e-8;

% betas=0.6655;
betas=0.6;
alpha=8*log((sqrt(1-betas)+1)/sqrt(betas));
ls=1/(alpha*sqrt(1-betas));
psim=log(betas)/alpha;

options2 = optimset('TolFun',tol,'Tolx',tol,'Display','off');

%%
close all;

fpsi=@(psi)1./sqrt(1.0-betas*exp(-alpha*psi));

dr=0.02; rmax=1.5/1;
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
fP=@(psi)betas*exp(-alpha*psi);
Pr=betas*exp(-alpha*psi);
fdP=@(psi)-betas*exp(-alpha*psi).*(alpha);
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
title(['\beta_s=',num2str(betas),', l_s=',num2str(ls)]);
subplot(222);
plot(rr,Pr,'linewidth',2);
xlabel('r');ylabel('P');
title(['\alpha=',num2str(alpha),', l_{sp}=',num2str(lsp,3)]);
subplot(223);
plot(rr,Jt,'linewidth',2);
xlabel('r');ylabel('J_\theta');
title(['h=',num2str(h)]);
subplot(224);
plot(rr,Bz,'linewidth',2);
xlabel('r');ylabel('B_z');
title(['\psi_m=',num2str(psim)]);

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
  'PaperSize',[screenposition(3:4)]);
print(gcf,'-dpng',['RR_betas=',num2str(betas),'_ls=',num2str(ls),...
    '_alpha=',num2str(alpha),'.png']);

datq_rr=[rr;Pr;psi;Bz;Jt].';
save('rr_datq.mat','datq_rr','betas','ls');