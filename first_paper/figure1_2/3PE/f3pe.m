% 2020-07-31 15:08, Hua-sheng Xie, ENN
% Equation for solve 3PE parameters, (betas,rs,Lp)->(alpha,sigma,psim)
function F=f3pe(x)
global betas ls alpha pm tol alpha0 ls0 psim0 psig;

sigma=x;

% pm=1;
alpha=((16*sigma*ls*sqrt(1-betas)+7)+pm*sqrt((16*sigma*ls*sqrt(1-betas)+7)^2-...
    4*sigma*48*ls*sqrt(1-betas)))/(2*48*ls*sqrt(1-betas));

% F0=@(alpha)ls*16*alpha*(3*alpha-sigma)*sqrt(1-betas)+(7*alpha-sigma);
% F0(alpha)


fP=@(psi)betas-8*betas*alpha^2/(7*alpha-sigma).*(...
    (psi<0).*(8*psi.*(1+4*sigma.*psi)+...
    (sigma+alpha).*((1-exp(16*alpha.*psi))/(8*alpha^2)))+...
    (psi>=0).*((1-exp(-8*alpha.*psi))/(alpha)-...
    (sigma+alpha).*(1-exp(-16*alpha.*psi))/(8*alpha^2) ));

F1=@(psim)fP(psim)-1+1e-6; %20-08-01 21:25
% F11=@(psim)(F1(psim)).^2;

options2 = optimset('TolFun',tol,'Tolx',tol,'Display','off');
psim=fsolve(F1,psig,options2);
% psim=fzero(F1,psig,options2);
% psim=fsolve(F1,psim0,options2);
% psim=fzero(F1,psim0,options2);
% psim=fminsearch(F11,psig,options2);

% psim=real(psim);

fpsi=@(psi)1./sqrt(1.0-fP(psi));

F=integral(fpsi,0,psim,'RelTol',tol,'AbsTol',tol)+1/4;

% F=real(F);  %%

end
