% 2020-07-28 21:05, Hua-sheng Xie, ENN
% Solve MRR parameters, (betas,ls)->(alpha,sigma,psim)
% 20-07-29 16:51 update
% 20-08-02 11:32 update from solve alpha to solve simga
close all;clear;clc;
format long;
mu0=4e-7*pi;

global betas ls n tol alpha0 ls0 psim0 psig;
n=2;
tol=1e-10;
icoil=1;

if (icoil==1) %five coils
    betas=0.7;
    Be=0.5;
    Rs=0.17;
    Rw=0.4;
    fprintf('icoil= %2d\n',icoil);
elseif (icoil==2)   %three coils
    betas=0.7;
    Be=1;
    Rs=0.075;
    Rw=0.14;
    fprintf('icoil= %2d\n',icoil);
end


alpha0=8*log((sqrt(1-betas)+1)/sqrt(betas));
ls0=1/(alpha0*sqrt(1-betas));
psim0=log(betas)/alpha0;
ls=0.5*ls0;

psig=-0.045;
options = optimset('TolFun',tol,'Tolx',tol,'Maxfuneval',10000,'MaxIter',5000);

x0=3.7;
x=fsolve(@fmrr0,x0,options)
F=fmrr0(x)
sigma=real(x);
alpha=((1/ls)/sqrt(1-betas)+n*sigma);
F1=@(psim)exp(-alpha*psim)*(sigma*psim+1)^n-1/betas;

options2 = optimset('TolFun',tol,'Tolx',tol,'Display','off');
psim=fsolve(F1,psig,options2);
F1=F1(psim);

%%
close all;
fpsi=@(psi)1./sqrt(1.0-betas*exp(-alpha*psi).*(sigma*psi+1).^n);

dr=0.005; rmax=Rw/Rs;
rr=0:dr:rmax;
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
close all;
figure('unit','normalized','DefaultAxesFontSize',16,...
    'DefaultAxesFontWeight','bold',...
    'DefaultAxesLineWidth',2,...
    'position',[0.5,0.5,0.5,0.6]);
fP=@(psi)betas*exp(-alpha*psi).*(sigma*psi+1).^n;
Pr=betas*exp(-alpha*psi).*(sigma*psi+1).^n;
fdP=@(psi)betas*exp(-alpha*psi).*(sigma*psi+1).^(n-1).*(-alpha*(sigma*psi+1)+n*sigma);
dP=betas*exp(-alpha*psi).*(sigma*psi+1).^(n-1).*(-alpha*(sigma*psi+1)+n*sigma);
Jt=-rr.*dP;
Bz=-sqrt(1-Pr).*(rr<=1/sqrt(2))+sqrt(1-Pr).*(rr>1/sqrt(2));

ids=find(abs(rr-1)==min(abs(rr-1)));ids=ids(1);
lsp=real(-(Pr(ids)/((Pr(ids+1)-Pr(ids-1))/(rr(ids+1)-rr(ids-1)))));

ind=find(rr<=1.0);
h=real(fdP(psim)/sum(2*dP(ind).*rr(ind)*dr));
ave_beta=sum(2*Pr(ind)./max(Pr).*rr(ind)*dr)/sum(2.*rr(ind)*dr);
tot_I=pi*sum(2*Jt(ind).*rr(ind)*dr);
tot_I2=pi*sum(2*Jt.*rr*dr);
Izeta0=tot_I2*((Be*Rs)/(2*mu0)); %include units:2020.10.23 15:51
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('************************** \n');
fprintf('Izeta0= %17.16f \n',Izeta0);


subplot(221);
plot(rr*Rs*100,psi.*(Be*Rs^2),'linewidth',2);hold all;
set(gca,'Fontsize',15);
xlabel('r');ylabel('\psi');xlim([0 100*Rw]);
title(['\beta_s=',num2str(betas),', l_s=',num2str(ls),', l_{sp}=',num2str(lsp,3)]);

subplot(222);
plot(rr*Rs*100,Pr.*(Be^2/(2*mu0)),'linewidth',2);
set(gca,'Fontsize',15);
xlabel('r');ylabel('P');xlim([0 100*Rw]);
title(['l_{s0}=',num2str(ls0),', \alpha=',num2str(alpha),', \psi_g=',num2str(psig)]);

subplot(223);
plot(rr*Rs*100,Jt.*(Be/(2*mu0*Rs)),'linewidth',2);
set(gca,'Fontsize',15);
xlabel('r');ylabel('J_\theta');xlim([0 100*Rw]);
title(['\sigma=',num2str(sigma),', h=',num2str(h),', n=',num2str(n)]);

subplot(224);
plot(rr*Rs*100,Bz.*Be,'linewidth',2);
set(gca,'Fontsize',15);
xlabel('r');ylabel('B_z');xlim([0 100*Rw]);
title(['\psi_m=',num2str(psim),', \alpha_0=',num2str(alpha0),', x0=',num2str(x0)]);
 
psibc =max(psi.*(Be*Rs^2));
fprintf('F1= %17.16f \n',F1);

C1s=alpha/(Be*Rs^2);
C2s=sigma/(Be*Rs^2);

fprintf('psibc= %17.16f \n',psibc);
psi1d=psi;  pr1d=Pr;    j1d=Jt;     Bz1d=Bz;    rr1d=rr;

if (1==1)
    save('MRR1d.mat','betas','ls','alpha','sigma','psibc','n','Rs','Rw','psi1d','pr1d',...
        'j1d','Bz1d','rr1d','C1s','C2s','n','Be','Rs','Rw','Izeta0','icoil');
end



