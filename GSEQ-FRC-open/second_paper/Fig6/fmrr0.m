% 2020-07-28 21:08, Hua-sheng Xie, ENN
% Equation for solve MRR parameters, (betas,rs,Lp)->(alpha,sigma,psim)
% 20-08-02 11:32 update from solve alpha to solve simga
function F=fmrr0(x)
% global betas ls n tol;
global betas ls n tol alpha0 ls0 psim0 psig;

sigma=x;

% n=2;
alpha=((1/ls)/sqrt(1-betas)+n*sigma);

F1=@(psim)exp(-alpha*psim)*(sigma*psim+1)^n-1/betas;
% F1=@(psim)betas*exp(-alpha*psim).*(sigma*psim+1).^n-1-1e-6;

options2 = optimset('TolFun',tol,'Tolx',tol,'Display','off');
psim=fsolve(F1,psig,options2);
% psim=fzero(F1,-0.005,options2);
% psim=fsolve(F1,psim0,options2);
% psim=fzero(F1,psim0,options2);

psim=real(psim);

fpsi=@(psi)1./sqrt(1.0-betas*exp(-alpha*psi).*(sigma*psi+1).^n);

F=integral(fpsi,0,psim,'RelTol',tol,'AbsTol',tol)+1/4;

% F=real(F);  %%

end
