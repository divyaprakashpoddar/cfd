function [aw, ap, ae, su] = get_coefficients(N, delX, P, gamma, area_cv, source, T)
    % Define the east and west grid points
    %W = P + 1;
    %E = P - 1;
    %delX_WP = myGrid(W) - myGrid(P);
    %delX_PE = myGrid(P) - myGrid(E);
    
    if(P == 1)
        aw = 0;
        ae = gamma*area_cv/delX;
        ap = aw + ae - (source(2) - 2*gamma*area_cv/delX);
        su = source(1) + 2*gamma*area_cv*T(1)/delX;
        return;
    elseif (P == N)
        aw = gamma*area_cv/delX;
        ae = 0;
        ap = aw + ae - (source(2) - 2*gamma*area_cv/delX);
        su = source(1) + 2*gamma*area_cv*T(end)/delX;
        return;
    end
    
    aw = gamma*area_cv/delX;
    ae = gamma*area_cv/delX;
    ap = aw + ae - source(2);
    su = source(1);
    
end