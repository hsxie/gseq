% 2020-07-31 10:58, Hua-sheng Xie, ENN
% Equation for solve MRR2 parameters, (betas,rs,Lp)->(alpha,sigma,psim)
% n>1
function F=fmrr2(x)
global betas ls alpha n tol alpha0 ls0 psim0 psig;

sigma=x;

F1=@(psim)exp(-alpha*psim-sigma*psim^n)-1/betas;

options2 = optimset('TolFun',tol,'Tolx',tol,'Display','off');
psim=fsolve(F1,psig,options2);
% psim=fzero(F1,-0.005,options2);
% psim=fsolve(F1,psim0,options2);
% psim=fzero(F1,psim0,options2);

psim=real(psim);

fpsi=@(psi)1./sqrt(1.0-betas*exp(-alpha*psi-sigma*psi.^n));

F=integral(fpsi,0,psim,'RelTol',tol,'AbsTol',tol)+1/4;

% F=real(F);  %%

end
