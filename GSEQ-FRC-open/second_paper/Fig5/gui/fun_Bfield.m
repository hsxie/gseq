% Hua-sheng Xie, CCF-ENN, huashengxie@gmail.com, 2019-06-16 15:34
function [psi,Br,Bz]=fun_Bfield(r0,z0,rl,zl,hl,wl,Il,Nl,imethod)
mu0=2.0e-7; % =mu0/(2*pi)
psi=0.0;
Br=0.0;
Bz=0.0;
nl=length(rl);

r0=r0+1e-8;
if(imethod==1) % Green function
    for jl=1:nl % to update
        k2=4*r0*rl(jl)/((r0+rl(jl))^2+(z0-zl(jl))^2);
        [EK2,EE2]=fun_ellipEK(k2,1);
        psi=psi+mu0*Nl(jl)*Il(jl)/(sqrt(k2))*sqrt(r0*rl(jl))*((2-k2)*EK2-2*EE2);
        Br=Br+mu0*Nl(jl)*Il(jl)*(z0-zl(jl))/r0...
            /sqrt((r0+rl(jl))^2+(z0-zl(jl))^2)*((r0^2+...
            rl(jl)^2+(z0-zl(jl))^2)/((r0-rl(jl))^2+(z0-zl(jl))^2)*EE2-EK2);
        Bz=Bz+mu0*Nl(jl)*Il(jl)...
            /sqrt((r0+rl(jl))^2+(z0-zl(jl))^2)*((rl(jl)^2-r0^2 ...
            -(z0-zl(jl))^2)/((r0-rl(jl))^2+(z0-zl(jl))^2)*EE2+EK2);
    end
    
elseif(imethod==2)
    
    for jl=1:nl % to update
        
        rmin=rl(jl)-0.5*hl(jl);rmax=rl(jl)+0.5*hl(jl);
        zmin=zl(jl)-0.5*wl(jl);zmax=zl(jl)+0.5*wl(jl);
        
        nr=2^4+1;nz=2^8+1;
        dr=(rmax-rmin)/(nr-1);
        dz=(zmax-zmin)/(nz-1);
        rr=rmin:dr:rmax;
        zz=zmin:dz:zmax;
        I1=Nl(jl)*Il(jl)/(nr*nz);
        for jr=1:nr
            for jz=1:nz
                r1=rr(jr);
                z1=zz(jz);
        
        k2=4*r0*r1/((r0+r1)^2+(z0-z1)^2);
        [EK2,EE2]=fun_ellipEK(k2,1);
        psi=psi+mu0*I1/(sqrt(k2))*sqrt(r0*r1)*((2-k2)*EK2-2*EE2);
        Br=Br+mu0*I1*(z0-z1)/r0...
            /sqrt((r0+r1)^2+(z0-z1)^2)*((r0^2+...
            r1^2+(z0-z1)^2)/((r0-r1)^2+(z0-z1)^2)*EE2-EK2);
        Bz=Bz+mu0*I1...
            /sqrt((r0+r1)^2+(z0-z1)^2)*((r1^2-r0^2 ...
            -(z0-z1)^2)/((r0-r1)^2+(z0-z1)^2)*EE2+EK2);
            end
        end
    end

else % Yuan Baoshan's book p110
    R1=r0; Z1=z0;
    for jl=1:nl % to update
    
    A1=rl(jl);
    B1=zl(jl);
    C1=hl(jl)*0.5;
    D1=wl(jl)*0.5;
    
    A=[1,48188,21101];
    C=[A1-C1,B1-D1,0.0];
    D=[A1+C1,B1+D1,0.5*pi];
    N=57091;
    X=0.*C;
    
    for M=1:N
        G=1.0;
        for L=1:3
            Z=M*A(L)/N;
            U=Z-floor(Z);
            V=(C(L)-D(L))*U;
            W=(U-1.0)*V;
            X(L)=C(L)+(2.0*W-V)*U;
            G=6.0*W*G;
        end
        E1=4.0*R1*X(1)/((R1+X(1))^2+(Z1-X(2))^2);
%         E2=sqrt(1.0-E1*sin(X(3)^2)); % wrong, 19-06-17 13:18
        E2=sqrt(1.0-E1*sin(X(3))^2);
        E3=(Z1-X(2))/(R1+1e-9)/sqrt((R1+X(1))^2+(Z1-X(2))^2);
        E4=(R1^2+X(1)^2+(Z1-X(2))^2)/((R1-X(1))^2+(Z1-X(2))^2);
        E5=2e-7*Nl(jl)*Il(jl)/(4*C1*D1);
        F1=E5*E3*(-1.0/E2+E4*E2);
        Br=Br+F1*G/N;
        E6=1.0/sqrt((R1+X(1))^2+(Z1-X(2))^2);
        E7=(R1^2-X(1)^2+(Z1-X(2))^2)/((R1-X(1))^2+(Z1-X(2))^2);
        F2=E5*E6*(1.0/E2-E7*E2);
        Bz=Bz+F2*G/N;
        
        % psi add by Huasheng Xie, 19-06-17 12:25
        E8=sqrt((R1+X(1))^2+(Z1-X(2))^2);
        E9=1.0-0.5*E1;
        F3=E5*E8*(E9/E2-E2);
        psi=psi+F3*G/N;
    end

    end

end
  
end