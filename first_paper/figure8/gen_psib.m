% CCF-ENN, huashengxie@gmail.com, 19-06-26 14:23
close all; clear; clc;
icase=8;

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
elseif(icase==7) % one frc, 19-10-07 09:06 for compress
zw=3.0;rw=0.7;rend=0.4;
elseif(icase==8) % one frc, 19-10-07 09:06 for compress
zw=0.75;rw=0.17;rend=0.15;
end

zz=-zw:0.005*zw:zw;
zb=zz;
rb=rw+0.0*zb;

% psiw=0.00376; psim=0.5*psiw;
% psim=0.00376; psiw=2.0*psim;
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
    zm1=-0.8*zw+zc; zm2=0.8*zw+zc;  psiw2=1.0*psiw;
    psib(ind1)=0.5*(psim-psiw)*(tanh(20*(zz(ind1)-zm1))+1.0)+psiw;
    psib(ind2)=0.5*(psim-psiw2)*(tanh(20*(zm2-zz(ind2)))+1.0)+psiw2;
elseif(icase==7)
    psim=0.00376; psiw=2.0*psim;
    zcon1=0.85; zcon2=zcon1+0.1;
    zc=(zcon1-0.0)*zw;   zc1=-(zcon1-0)*zw;
    zm1=-zcon2*zw; zm2=zcon2*zw;

    ind1=find(abs(zz)<=zc);
    ind2=find(abs(zz)>=zm2);
    ind3=find(zz>=zc & zz<zm2); %right
    ind4=find(zz>=zm1 & zz<zc1); %right 
    rb(ind1)=1.0*rw;    psib(ind1)=psim;
    rb(ind2)=1.0*rend;  psib(ind2)=psim;
    rb(ind3)=rw+(rend-rw)/(zm2-zc)*(zz(ind3)-zc);       psib(ind3)=psim;
    rb(ind4)=rend+(rw-rend)/(zc1-zm1)*(zz(ind4)-zm1);   psib(ind4)=psim;
elseif(icase==8)
    psim=1.0*0.006384168;
%     psim=0.0091591;
%     zcon1=0.6; zcon2=zcon1+0.05/1.5;
    zcon1=0.6; zcon2=zcon1+0.05/0.75;
    zc=(zcon1-0.0)*zw;   zc1=-(zcon1-0)*zw;
    zm1=-zcon2*zw; zm2=zcon2*zw;

    ind1=find(abs(zz)<=zc);
    ind2=find(abs(zz)>=zm2);
    ind3=find(zz>=zc & zz<zm2); %right
    ind4=find(zz>=zm1 & zz<zc1); %right 
    rb(ind1)=1.0*rw;    psib(ind1)=psim;
    rb(ind2)=1.0*rend;  psib(ind2)=psim;
    rb(ind3)=rw+(rend-rw)/(zm2-zc)*(zz(ind3)-zc);       psib(ind3)=psim;
    rb(ind4)=rend+(rw-rend)/(zc1-zm1)*(zz(ind4)-zm1);   psib(ind4)=psim;    
    
elseif(icase==4)
    psim=-0.00376; psiw=2.0*psim;
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
    
%     psib(ind1)=0.5*(psim-psiw)*(tanh(20*(zz(ind1)-zm1))+1.0)+psiw;
%     psib(ind2)=0.5*(psim-psiw2)*(tanh(20*(zm2-zz(ind2)))+1.0)+psiw2;
%     psib(ind3)=0.5*(psim-psiw)*(tanh(20*(zz(ind3)-zm3))+1.0)+psiw;
    
else
    psib=psim*(1+cos(3*pi*zz/zw)+(1.0*zz./zw).^2)+psiw;
end

hf=figure('unit','normalized','Position',[0.6 0.3 0.4 0.2],...
   'DefaultAxesFontSize',15); 
if icase==7 || icase==8
    plot(zz,rb,'.-m');grid on;
    set(gca,'Fontsize',23);
    xlabel('$z / Z_{max}$','interpreter','latex','FontSize',17);
    ylabel('$r / R_{w}$','interpreter','latex','FontSize',17);
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
    h=legend('$R_{w}(z)$');set(h,'interpreter','latex');
    legend('boxoff');
    ylim([0,rw])
else 
    plot(zb,psib,'linewidth',2);grid on;
    set(gca,'Fontsize',17);
    xlim([-zw,zw]);
    ylim([0,2*psiw]);
    xlabel('$Z \ (m)$','interpreter','latex','FontSize',17);
    ylabel('$\psi \ (Wb)$','interpreter','latex','FontSize',17);
    set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);

end


save('psib.mat','zb','rb','psib');

% h=legend('$case \ 1$');
% set(h,'interpreter','latex');
