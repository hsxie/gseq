% CCF-ENN, huashengxie@gmail.com, 19-06-26 14:23
close all; clear; clc;
icase=7;

figure('position',[50 100 1000 400]);
set(gca,'nextplot','replacechildren');
set(gcf,'DefaultAxesFontSize',15);

if(icase==1)
    zw=1.0;rw=0.17;
elseif(icase==2)
    zw=1.5;rw=0.17;
elseif(icase==3) % one frc, for acceleration and translation 
    zw=1.0;rw=0.17;
elseif(icase==4)
    zw=1.0;rw=0.17;
elseif(icase==5) % two frc, 19-10-04 20:24, for collide
    zw=1.5;rw=0.17;
elseif(icase==6) % one frc, 19-10-07 09:06 for compress
    zw=1.0;rw=0.17;
elseif(icase==7) % one frc, grass_ft psi_bound 
    zw=1.5;rw=0.17; 
    psi_m1 =-4;psi_wall =-1;psi_m2 =-4;
    Z_m1 =-0.5;Z_m2 =0.5;Z_c1 =-0.25;Z_c2=0.25;  
    
    mo =385;m=129;
    rm1 =4.0;rm2 =4.0;
    jz1 =129+((m-1)/4);
    jz2 =257-((m-1)/4);
elseif(icase==8) %
    zw=1.5;rw=1.0;    
end

% zz=-zw:0.005*zw:zw;
zz =linspace(-zw,zw,385);
zb=zz;
rb=rw+0.0*zb;

if(icase==2)
    psim=0.00376; psiw=4.0*psim;
    zc=-0.1*zw;

    ind1=find(zz<zc);
    ind2=find(zz>=zc);
    zm1=-0.3*zw+zc; zm2=0.3*zw+zc;  psiw2=1.0*psiw;
    % psib(ind1)=0.5*(psiw-psim)*(tanh(2*(zz(ind1)+zm))+1.0)+psim;
    % psib(ind2)=0.5*(psiw-psim)*(tanh(2*(zm-zz(ind2)))+1.0)+psim;
    psib(ind1)=0.5*(psim-psiw)*(tanh(20*(zz(ind1)-zm1))+1.0)+psiw;
    psib(ind2)=0.5*(psim-psiw2)*(tanh(20*(zm2-zz(ind2)))+1.0)+psiw2;
elseif(icase==1)
    psim=0.00376; psiw=2.0*psim;
    zc=-0.4*zw;

    ind1=find(zz<zc);
    ind2=find(zz>=zc);
    zm1=-0.2*zw+zc; zm2=0.2*zw+zc;  psiw2=1.0*psiw;
    psib(ind1)=0.5*(psim-psiw)*(tanh(20*(zz(ind1)-zm1))+1.0)+psiw;
    psib(ind2)=0.5*(psim-psiw2)*(tanh(20*(zm2-zz(ind2)))+1.0)+psiw2;
elseif(icase==3)
    psim=0.00376; psiw=2.0*psim;
    zc=-0.2*zw;

    ind1=find(zz<zc);
    ind2=find(zz>=zc);
    zm1=-0.2*zw+zc; zm2=0.2*zw+zc;  psiw2=1.0*psiw;
    psib(ind1)=0.5*(psim-psiw)*(tanh(20*(zz(ind1)-zm1))+1.0)+psiw;
    psib(ind2)=0.5*(psim-psiw2)*(tanh(20*(zm2-zz(ind2)))+1.0)+psiw2;
elseif(icase==6)
    psim=0.00376; psiw=2.0*psim;
    zc=-0.0*zw;

    ind1=find(zz<zc);
    ind2=find(zz>=zc);
    zm1=-0.3*zw+zc; zm2=0.3*zw+zc;  psiw2=1.0*psiw;
    psib(ind1)=0.5*(psim-psiw)*(tanh(20*(zz(ind1)-zm1))+1.0)+psiw;
    psib(ind2)=0.5*(psim-psiw2)*(tanh(20*(zm2-zz(ind2)))+1.0)+psiw2;
%     psib(ind1)=-0.5*(psim-psiw)*(tanh(20*(zz(ind1)-zm1))+1.0)+psiw;
%     psib(ind2)=-0.5*(psim-psiw2)*(tanh(20*(zm2-zz(ind2)))+1.0)+psiw2;
% psib(ind2)=0.;
elseif(icase==4)
    psim=0.00376; psiw=2.0*psim;
    zc=-0.0*zw;

    ind1=find(zz<zc);
    ind2=find(zz>=zc);
    zm1=-0.2*zw+zc; zm2=0.2*zw+zc;  psiw2=1.0*psiw;
    psib(ind1)=0.5*(psim-psiw)*(tanh(20*(zz(ind1)-zm1))+1.0)+psiw;
    psib(ind2)=0.5*(psim-psiw2)*(tanh(20*(zm2-zz(ind2)))+1.0)+psiw2;
elseif(icase==5) % two frc, 19-10-04 20:24
    psim=0.00376; psiw=2.0*psim;
    zc=-0.4*zw;
    zm1=-0.15*zw+zc; zm2=0.15*zw+zc; zm3=-0.15*zw-zc; zm4=0.15*zw-zc;  
    psiw2=1.0*psiw;
    psib=0.*zz;
    for j=1:length(zz)
        if(zz(j)<zc)
            psib(j)=0.5*(psim-psiw)*(tanh(20*(zz(j)-zm1))+1.0)+psiw;
        elseif(zz(j)>=zc && zz(j)<0)
            psib(j)=0.5*(psim-psiw2)*(tanh(20*(zm2-zz(j)))+1.0)+psiw2;
        elseif(zz(j)>=0 && zz(j)<-zc)
            psib(j)=0.5*(psim-psiw2)*(-tanh(20*(zm3-zz(j)))+1.0)+psiw2;
        else
            psib(j)=0.5*(psim-psiw)*(-tanh(20*(zz(j)-zm4))+1.0)+psiw;
        end   
     
    end
elseif(icase==7) % one frc, grass_ft psi_bound   
    sub_plot=2;
    if sub_plot ==1
        for i=1:size(zz,2)
           if zz(i)<=Z_m1
               psib(i) =psi_m1;
           elseif zz(i)>Z_m1 && zz(i)<Z_c1 
               psib(i) =0.5*(psi_m1+psi_wall)+0.5*(psi_m1-psi_wall)*cos(pi*(zz(i)-Z_m1)/(Z_m1-Z_c1));
           elseif zz(i)<Z_c2 && zz(i)>=Z_c1    
                psib(i) =psi_wall;
           elseif zz(i)<Z_m2 && zz(i)>=Z_c2  
               psib(i) = 0.5*(psi_m2+psi_wall)+0.5*(psi_wall-psi_m2)*cos(pi*(zz(i)-Z_c2)/(Z_m2-Z_c2));
           elseif zz(i)>=Z_m2
               psib(i) =psi_m2;
           end    
        end
    elseif sub_plot ==2
        for i=1:mo
            if (i<=jz1 && i>129) 	%two
                psib(i) =(-0.5*(rm1+1.)-0.5*(rm1-1.)*cos(pi/(zz(jz1)+0.5)*(zz(i)+0.5)));
            elseif (i>jz1 && i<=jz2)   %three
                psib(i) =-1.0;
            elseif(i>jz2 && i<=257) 	%four
                psib(i) =(-0.5*(rm2+1.)-0.5*(rm2-1.)*cos(pi/(zz(jz2)-0.5)*(zz(i)-0.5)));
            elseif (i<=129)    %one
                psib(i) =(-0.5*(rm1+1.)-0.5*(rm1-1.)*cos(pi/(zz(jz1)+0.5)*(zz(129)+0.5)));
            elseif (i>=257)    %five
                psib(i) =(-0.5*(rm2+1.)-0.5*(rm2-1.)*cos(pi/(zz(jz2)-0.5)*(zz(257)-0.5)));
            end
        end
    end

elseif(icase==8)
    for i=1:size(zz,2)
        psib(i) =0.0;
    end
else
    psib=psim*(1+cos(3*pi*zz/zw)+(1.0*zz./zw).^2)+psiw;
end

plot(zb,psib,'m','linewidth',2);
xlim([-zw,zw]);
grid on;box on;
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',0.5);
xlabel('$Z(m)$','interpreter','latex','FontSize',17);
ylabel('$\psi$','interpreter','latex','FontSize',17);


save('psib.mat','zb','rb','psib');
