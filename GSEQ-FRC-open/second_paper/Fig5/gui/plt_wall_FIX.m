% 19-06-20 07:50, plot wall
zw=3.5;rw=1;

zz=-zw:0.01*zw:zw;

ind1=find(zz<0);
ind2=find(zz>=0);
% a=0.15; zm=0.8*zw; ra=0.4*rw; rb=0.8*rw;
% a=0.15; zm=0.75*zw; ra=0.4*rw; rb=0.7*rw;
a=0.15; zm=0.6*zw; ra=0.3*rw; rb=0.4*rw;
% rr(ind1)=0.15*(tanh(5*(zz(ind1)+0.8*zw))+1.0)+0.4;
wall=1;
if(wall==1)
    rr(ind1)=0.5*(rb-ra)*(tanh(5*(zz(ind1)+zm))+1.0)+ra;
    rr(ind2)=0.5*(rb-ra)*(tanh(5*(zm-zz(ind2)))+1.0)+ra;
else
    rr=0.15+0.*zz;
end
plt=2;
if(plt==1)
    figure;
    subplot(211);
    xlabel('z');ylabel('r'); hold on;
else
    set(gcf,'currentaxes',hAxes1);
end

plot(zz,rr,'r--','linewidth',2); hold on;
plot(zz,-rr,'r--','linewidth',2);

if(plt==1)
    subplot(212);
    xlabel('z');ylabel('\psi');
elseif(plt==2)
    set(gcf,'currentaxes',hAxes2);
    cla;
    xlabel('z');ylabel('\psi at wall'); hold on;
    grid on;
    xlim([-zw,zw]);
end

nz=length(zz);
psi=0.*zz;
Br=0.*zz;
Bz=0.*zz;
for jz=1:nz
    [psi(jz),Br(jz),Bz(jz)]=fun_Bfield(rr(jz),zz(jz),rl,zl,hl,wl,Il,Nl,imethod);
end

if(plt~=3)
    hold on;
    plot(zz,psi,'linewidth',2);
end

rb=rr;zb=zz;psib=psi;
save('psib_FIX.mat','zb','rb','psib');
save('../psib_FIX.mat','zb','rb','psib','Bz');
Z0=(size(psi,2)+1)/2;
% fprintf('psib=%17.16f \n',psi(65));

No=find(zz>-3 & zz<0);
[Bz_min_no]=find(Bz==min(Bz(No)));
fprintf('the position of mid-plane of two FRC: %5.4d \n',Bz_min_no);
fprintf('the axial position of mid-plane of two FRC: %5.4d \n',zz(Bz_min_no));
fprintf('the magnetic field of mid-plane of two FRC: %5.4d \n',Bz(Bz_min_no));


jN0=floor(0.5*nz)-floor(1.5/(0.01*zw));
fprintf('psib=%17.16f \n',psi(jN0));


