% Hua-sheng Xie, CCF-ENN, huashengxie@gmail.com, 2019-06-06 16:34
% Elliptic functions for G-S equilibrium solver's Green function of coils
% [E,K]=fun_ellipEK(k2,apx) % 19-06-07 11:05 [E,K] should be [K,E]!
function [K,E]=fun_ellipEK(k2,apx)
%   apx=1;
  if(nargin==1)
    apx=1;
  end
  
  if(apx==2)
    [K,E]=ellipke(k2,1e-8);
  elseif(apx==0)
    K=ellipticK(k2);
    E=ellipticE(k2);
  elseif(apx==1)
    % from http://people.sc.fsu.edu/~jburkardt%20/f_src/special_functions/special_functions.f90
    % Accurate and fast
    K=0.0.*zeros(length(k2),1);
    E=0.0*K;
    for j=1:length(k2)
      pk = 1.0D+00 - k2(j);

      if ( k2 == 1.0 )
        K(j) = 1.0D+300;
        E(j) = 1.0D+00;
      else

        ak = (((0.01451196212D+00   * pk + 0.03742563713D+00 ) * pk + 0.03590092383D+00 ) * pk + ...
            0.09666344259D+00 ) * pk  + 1.38629436112D+00;

        bk = (((  0.00441787012D+00   * pk+ 0.03328355346D+00 ) * pk  + 0.06880248576D+00 ) * pk ...
          + 0.12498593597D+00 ) * pk  + 0.5D+00;

        K(j) = ak - bk * log ( pk );

        ae = ((( 0.01736506451D+00   * pk + 0.04757383546D+00 ) * pk  + 0.0626060122D+00  ) * pk ...
          + 0.44325141463D+00 ) * pk + 1.0D+00;

        be = (((  0.00526449639D+00   * pk  + 0.04069697526D+00 ) * pk  + 0.09200180037D+00 ) * pk ...
          + 0.2499836831D+00  ) * pk;

        E(j) = ae - be * log ( pk );

      end
    end
  
  else % Taylor expansion for Elliptic function
    K=0.5*pi*(1+k2/4+9/64*k2.^2+(15/48)^2*k2.^3);
    E=0.5*pi*(1-k2/4-3/64*k2.^2-45/48^2*k2.^3); 
  end
  
end