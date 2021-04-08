% Hua-sheng Xie, CCF-ENN, huashengxie@gmail.com, 2019-06-05 15:36
% Free boundary G-S equilibrium for FRC with point coils
% Green function method
% 19-06-07 14:25 seems ok, but slow and senstive to profiles
% 20-09-01 12:51 test MRR model, seems can run
close all;clear; clc;
irun=0; % =0, initial run; =1, from previous run
if(irun==0)
    clear; clc; irun=0;
    mu0=4*pi*1e-7;
    
    init=0;
    icase=2;
    jP=1; %jP=3; % choose P(psi) profile type
    rw=0.14; zw=0.75;
    
    if jP==1
        load('MRR1d.mat');
        Izeta0=-Izeta0;
        Pm0=Be^2/(2*mu0);
    end
    
    if icase==1
        hl=[0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05]*rw; % r width of coils
        rl=[0.7,0.8,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,1.1,0.8,0.7]*rw; % r position of coils
        zl=[-1.15,-1.1,-0.9,-0.6,-0.3,-0.1,0.0,0.1,0.3,0.6,0.9,1.1,1.15]*zw; % z position of coils
        wl=[0.1,0.1,0.1,0.3,0.3,0.1,0.1,0.1,0.3,0.3,0.1,0.1,0.1]*zw; % z width of coils
        Il=[10e5,5e5,1.0e5,1.0e5,1.0e5,1.0e5,0e4,1.0e5,1.0e5,1.0e5,1e5,5e5,10e5]*0.4; 
        % current of coils, Ampere
    elseif icase==2
        hl=[0.006,0.006,0.006];
        rl=[0.15,0.15,0.15];
        zl=[-0.7,0,0.7];
        wl=[0.1,1.2,0.1];
        Il=[1000,1000,1000]*1000;
    end
    
    nl=length(rl); % # of coils
    nR=2^7+1; nZ=2^8+1; % grid numbers
    
    Rlow=0; Rup=rw; Zlow=-1.0*zw; Zup=1.0*zw;
    rr=linspace(Rlow,Rup,nR); zz=linspace(Zlow,Zup,nZ);
    dR=rr(2)-rr(1); dZ=zz(2)-zz(1);
    [R,Z]=ndgrid(rr,zz); dR2=dR*dR; dZ2=dZ*dZ;
    Z0=(size(R,2)+1)/2;
    
%     Green function psi on the boundary
    nb=2*(nR-1)+(nZ-2);
    psib=zeros(nb,1); jrb=ones(nb,1); jzb=ones(nb,1);psiw_coil=zeros(nb,1);
    Br_coil=zeros(nb,1);Bz_coil=zeros(nb,1);
    
    for j=1:nb
        if(j<nR) % z=zmin
            jrb(j)=j+1; jzb(j)=1;
        elseif(j>(nb-nR+1)) % z=zmax
            jrb(j)=nb-j+2; jzb(j)=nZ;
        else % r=rmax
            jrb(j)=nR; jzb(j)=j-nR+2;
        end
        
        
        for jl=1:nl
            rmin=rl(jl)-0.5*hl(jl);rmax=rl(jl)+0.5*hl(jl);
            zmin=zl(jl)-0.5*wl(jl);zmax=zl(jl)+0.5*wl(jl);
            nr=10;nz=100;
            if (nr==1 || nz==1)
                rr_coil=rl(jl);
                zz_coil=zl(jl);
            else
                dr=(rmax-rmin)/(nr-1);
                dz=(zmax-zmin)/(nz-1);
                rr_coil=rmin:dr:rmax;
                zz_coil=zmin:dz:zmax;
            end

            I1=Il(jl)/(nr*nz);
            for jr=1:nr
                for jz=1:nz
                    r1=rr_coil(jr);
                    z1=zz_coil(jz);
                    k2=4*R(jrb(j),jzb(j))*r1/((R(jrb(j),jzb(j))+...
                        r1)^2+(Z(jrb(j),jzb(j))-z1)^2);
                    [EK2,EE2]=fun_ellipEK(k2,1);
                    if(k2~=0)
                        Grz_coil=1/(2*pi*sqrt(k2))*sqrt(R(jrb(j),...
                            jzb(j))*r1)*((2-k2)*EK2-2*EE2);
                        psiw_coil(j)=psiw_coil(j)+Grz_coil*mu0*I1;
                        
                    Br_coil(j)=Br_coil(j)+mu0*I1*(Z(jrb(j),jzb(j))-z1)/R(jrb(j),jzb(j))...
                        /sqrt((R(jrb(j),jzb(j))+r1)^2+(Z(jrb(j),jzb(j))-z1)^2)*((R(jrb(j),jzb(j))^2+...
                        r1^2+(Z(jrb(j),jzb(j))-z1)^2)/((R(jrb(j),jzb(j))-r1)^2+(Z(jrb(j),jzb(j))-z1)^2)*EE2-EK2);
                    Bz_coil(j)=Bz_coil(j)+mu0*I1/sqrt((R(jrb(j),jzb(j))+r1)^2+(Z(jrb(j),jzb(j))-z1)^2)*...
                        ((r1^2-R(jrb(j),jzb(j))^2-(Z(jrb(j),jzb(j))-z1)^2)/...
                        ((R(jrb(j),jzb(j))-r1)^2+(Z(jrb(j),jzb(j))-z1)^2)*EE2+EK2)/2/pi;
                    end
                end
            end
        end
%         psib(j)=psiw_coil(jrb(j),jzb(j));
        psib(j)=psiw_coil(j);
    end
    psib_coil=psib; % save for further use
    
    GrzJz=zeros(nR,nZ,nb);
    for jr=2:(nR-1)
        for jz=2:(nZ-1)
            for j=1:nb
                k2=4*R(jrb(j),jzb(j))*R(jr,jz)/((R(jrb(j),jzb(j))+...
                    R(jr,jz))^2+(Z(jrb(j),jzb(j))-Z(jr,jz))^2);
                [EK2,EE2]=fun_ellipEK(k2,1);
                if(k2~=0)
                    GrzJz(jr,jz,j)=1/(2*pi*sqrt(k2))*sqrt(R(jrb(j),...
                        jzb(j))*R(jr,jz))*((2-k2)*EK2-2*EE2);
                end
            end
        end
    end
    
    maxpsiw=max(abs(psib_coil));
    
    if(init==0)
        if jP==1
            C=1; C3=-0.0005;
            psi=0.01*(0.1-exp(-((R/rw-0.2).^2/0.3^2+(Z/zw).^2/0.4^2)))*maxpsiw;
            
        elseif jP==2
            C=0.0061; C3=-0.0005;
            psi=(0.5-exp(-((R/rw-0.3).^2/0.3^2+(Z/zw).^2/0.3^2)))*maxpsiw; % initial psi, to update
        end
        
        CI=1.0;
        CItmp=CI;
    else
        load('psiq.mat'); % 18-11-13 09:16
        psi=interp2(Zq,Rq,psiq,Z,R);
        CI=CIq;
    end
    Ctmp=C; Ctmp3=C3;
    psi0=psi; psi1=psi; psi2=psi;
    dpsim=max(max(abs(psi-psi0)))+0.1;
    
    it=1; % iterative times
    ita=0;
    itout=1;
    runtime0=0;
    
    psiout=psi;
    
else
    ita=itout;
    runtime0=runtime;
end
runtime=cputime;
nt=20000; % max iterative steps
ntout=30; % max out iterative steps
tol=1e-9; % tolerance for stop
dpsimout=100*tol;

%%
hf=figure('unit','normalized','Position',[0.5 0.1 0.5 0.5],...
    'DefaultAxesFontSize',15);

Be1=1.1*Be;
itCI=1;
wx=0.95;

% first loop
while(abs(Be1/Be-1)>0.0005 && itCI<200)
    itout=1;
    it=1;
while(itout<=ntout+ita && dpsimout>tol*2) % iteration, G-S eq
    itout
    
    Jzeta=0.*R;
    for jR=2:(nR-1)
        for jZ=2:(nZ-1)
            Jzeta(jR,jZ)=-1/R(jR,jZ)*((psi(jR+1,jZ)-2*psi(jR,jZ)+psi(jR-1,jZ))/(dR^2)...
                -1/R(jR,jZ)*(psi(jR+1,jZ)-psi(jR-1,jZ))/(2*dR)...
                +(psi(jR,jZ+1)-2*psi(jR,jZ)+psi(jR,jZ-1))/(dZ^2));
        end
    end
    Jzeta=Jzeta/mu0;
    Jztotal=sum(sum(Jzeta))*dR*dZ;
    
    if (1==0)
        psib=psib_coil;
        for jr=2:(nR-1) % update the boundary flux psi
            for jz=2:(nZ-1)
                for j=1:nb
                    psib(j)=psib(j)+mu0*CI*GrzJz(jr,jz,j)*Jzeta(jr,jz)*dR*dZ; 
                end
            end
        end
    else
        psib=CI*psib_coil;
        for jr=2:(nR-1) % update the boundary flux psi
            for jz=2:(nZ-1)
                for j=1:nb
                    psib(j)=psib(j)+mu0*CI*GrzJz(jr,jz,j)*Jzeta(jr,jz)*dR*dZ; 
                end
            end
        end
    end
    
    dpsim=10*tol; it=1;
    if(itout==1)
        jit=1; w=0.99; % parameter for iterative, =0.9 for default
        w0=w;
    elseif(itout>5)
%               w=0.85; w0=w;
    end
    
    %update psi loop
    while(it<=nt && dpsim>tol)     
        if(jP==1)
            dP=C*Pm0*betas.*exp(-alpha.*psi./(Be*Rs^2)).*(sigma.*psi./(Be*Rs^2)+1).^(n-1).*...
                (-alpha/(Be*Rs^2).*(sigma.*psi./(Be*Rs^2)+1)+n*sigma/(Be*Rs^2));
        elseif(jP==2) % for sech pressure profile, 19-06-07 11:20, suggested by DengBH
            C2=500;
            dP=-C2*C*tanh(C2*(psi-C3)).*sech(C2*(psi-C3))/mu0;
        end
        
        Jzeta=R.*dP;
        Jztotal=sum(sum(Jzeta))*dR*dZ;
        Jztotal_mid=0.0;
        
        for jz=Z0 % 19-06-23 22:00 update
            for jr=1:(jrb(jz)-1)
                Jztotal_mid=Jztotal_mid+2*pi*Jzeta(jr,jz)*R(jr,jz)*dR;
            end
        end
        
        P=C*Pm0*betas.*exp(-alpha.*psi./(Be*Rs^2)).*(sigma.*psi./(Be*Rs^2)+1).^n;
        
        if(itout>1)
            wc=w;
            if jP==1
                C=wc*C+(1-wc)*Ctmp*Izeta0/Jztotal_mid;
%                 C=wc*C+(1-wc)*Ctmp*Izeta/Jztotal;
            elseif jP==2
                C=wc*C+(1-wc)*Ctmp*Izeta/Jztotal;
            end
        end
        
        for i=2:nR-1
            for j=2:nZ-1
                if(jit==1) % one step iterative
                    psi(i,j)=w*psi0(i,j)+0.5*dZ2*( rr(i)/dR2*((psi0(i+1,j)-...
                        psi0(i,j))/(rr(i)+0.5*dR)-...
                        (psi0(i,j)-psi0(i-1,j))/(rr(i)-0.5*dR))+...
                        (psi0(i,j+1)+psi0(i,j-1))/dZ2+...
                        mu0*rr(i)^2*dP(i,j) )*(1-w); % 18-11-04 16:25, psi should be psi0
                elseif (jit==2)  % 3-point
                    psi(i,j)=1/(2/dR2+2/dZ2)*(...
                        (psi0(i+1,j)+psi0(i-1,j))/dR2-...
                        (psi0(i+1,j)-psi0(i-1,j))/(2*dR*rr(i))+...
                        (psi0(i,j+1)+psi0(i,j-1))/dZ2+...
                        mu0*rr(i)^2*dP(i,j));
                elseif (jit==3) % 5-point + 3-point
                    if((j>2 && j<nZ-1) && (i>2 && i<jrb(j)-1))
                        psi(i,j)=1/(30/(12*dR2)+30/(12*dZ2))*(...
                            (-psi0(i+2,j)+16*psi0(i+1,j)+16*psi0(i-1,j)-psi0(i-2,j))/(12*dR2)-...
                            (-psi0(i+2,j)+8*psi0(i+1,j)-8*psi0(i-1,j)+psi0(i-2,j))/(12*dR*rr(i))+...
                            (-psi0(i,j+2)+16*psi0(i,j+1)+16*psi0(i,j-1)-psi0(i,j-2))/(12*dZ2)+...
                            mu0*rr(i)^2*dP(i,j));
                    else
                        psi(i,j)=1/(2/dR2+2/dZ2)*(...
                            (psi0(i+1,j)+psi0(i-1,j))/dR2-...
                            (psi0(i+1,j)-psi0(i-1,j))/(2*dR*rr(i))+...
                            (psi0(i,j+1)+psi0(i,j-1))/dZ2+...
                            mu0*rr(i)^2*dP(i,j));
                    end
                end
            end
        end
        % boundary condition
        psi(1,:)=0.0;
        for j=1:nb
            psi(jrb(j),jzb(j))=psib(j);
        end
        
        it=it+1;
        dpsim=max(max(abs(psi-psi0)));
        psi2=psi1; psi1=psi0; psi0=psi;
        
        Ctmp=C;
    end
    
    dpsimout=max(max(abs(psi-psiout)));
    psiout=psi;
    
    if(mod(itout,1)==0) % plot
        fprintf('C = %5.4d \n',C);
        ido=find(psi==min(min(psi)));
        
        subplot(2,3,1:3); hold off;
        indtmp= psi>0; psio=psi; psio(indtmp)=NaN;
        indtmp=find(psi<0); psii=psi; psii(indtmp)=NaN;
        contour(Z,R,psio,5,'linewidth',2); hold on;
        contour(Z,R,psii,20,'linewidth',2); hold on;
        [ch,hh]=contour(Z,R,psi,[0,0],'r--','linewidth',3); hold on; % highligth separatix
        plot(Z(ido(1)),R(ido(1)),'mx','linewidth',3); % O-point
        text(Z(ido(1))+0.01*zw,R(ido(1)),'O point','fontsize',13);
        title(['(a) \Psi(R,Z)',...
            ', it=',num2str(itout),'/',num2str(nt),...
            ', max(\Psi^{n}-\Psi^{n-1})=',num2str(dpsim),...
            ', max(\Psi^{n}-\Psi^{n-1})_{out}=',num2str(dpsimout)]);
        
        ylabel(['R, r_w=',num2str(rw)]);
        xlabel(['I=',num2str(Jztotal,3),'Be=',num2str(Be,3),', Be1=',num2str(Be1,3),...
            ', CI=',num2str(CI,3),', itCI=',num2str(itCI)]);
        if(jP==1)
            strp='p''(\psi)=C(1-tanh(C2*psi-C3))';
        elseif(jP==2)
            strp='p(\psi)=C*sech(C2*(psi-C3))';
        end
        box on; %axis equal; 
        hold on;plot(zl,rl,'ro');
        for jl=1:nl
            text(zl(jl),rl(jl),['I_l=',num2str(Il(jl)/1e6),'MA']);
            hold on;
            if(Il(jl)>0)
                rectangle('position',[zl(jl)-0.5*wl(jl),rl(jl)-0.5*hl(jl),wl(jl),hl(jl)],...
                    'linewidth',1);hold on;                
            else
                rectangle('position',[zl(jl)-0.5*wl(jl),rl(jl)-0.5*hl(jl),wl(jl),hl(jl)],...
                    'linewidth',1);hold on;  
            end
        end
        rectangle('position',[Zlow,Rlow,(Zup-Zlow),(Rup-Rlow)]);
        hold off;
        
        subplot(234);
        psi_x=psi(:,Z0);
        itmp=floor(0.2*nR);
        rs=rr(find(abs(psi_x(itmp:end))==min(abs(psi_x(itmp:end))))+itmp-1);
        plot(rr,-Jzeta(:,Z0)/max(abs(Jzeta(:,Z0))),'m','linewidth',3);box on;grid on;
        set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
        xlabel('r (m)'); ylabel('J_{\theta}');
        title(['(b) J_m=',num2str(max(abs(Jzeta(:,Z0))),3),'; rs=',num2str(rs(1))]);
        
        subplot(235);
        E=0.5*(max(ch(1,:))-min(ch(1,:)))/rs(1);
        contourf(R,Z,dP); title(['(c) dP, tol=',num2str(tol),', E= ',...
            num2str(E)]);
        xlabel('R'); ylabel('Z'); box on;

        
        subplot(236);
        plot(rr,psi(:,Z0),'m','linewidth',3);box on;grid on;
        title(['(d) \Psi, w=',num2str(w),'; C=',num2str(C)]);
        xlabel('r'); ylabel('\psi');
        set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
        
        pause(0.0001);
    end
    itout=itout+1;
    
end

Be1=sqrt(2*mu0*max(P(:,Z0)));
% CI=CI*Be1/(0.05*Be1+0.95*Be);
%CI=CI*Be1/(0.05*Be1+0.95*Be) %SOR iteration
CI=wx*CI+(1-wx)*CItmp*Be/Be1;
itCI=itCI+1
CItmp=CI;
end
runtime=runtime0+cputime-runtime;

%%
Rq=R;Zq=Z;psiq=psi; CIq=CI;
save('psiq.mat','Rq','Zq','psiq','C','C3','CIq');
save('psiq1.mat','Rq','Zq','psiq','Jzeta','P','jrb','Z0','zl','wl',...
     'rl','hl','C','E','Bz_coil','psiw_coil','CI');

