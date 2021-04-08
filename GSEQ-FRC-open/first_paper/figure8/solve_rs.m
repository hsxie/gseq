function [rs]=solve_rs(psi,r)
    ind_min=find(psi==min(psi));
    rs=interp1(psi(ind_min(1):end-1),r(ind_min(1):end-1),0,'linear');
end