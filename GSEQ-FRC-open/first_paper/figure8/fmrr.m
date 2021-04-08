function F=fmrr(x)
global betas rs Lp;

alpha=x(1);
sigma=x(2);
psim=x(3);

n=2;

F=0.*x;
F(1)=(n*sigma-alpha)*sqrt(1-betas)+rs/Lp;
F(2)=exp(-alpha*psim)*(sigma*psim+1)^n-1/betas;

fpsi=@(psi)1./sqrt(1-betas*exp(-alpha*psi).*(sigma*psi+1).^n);
F(3)=-integral(fpsi,0,psim+1e-10)-1/4;


end