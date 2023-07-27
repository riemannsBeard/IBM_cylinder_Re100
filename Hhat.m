function [H, beta] = Hhat(grid, ib, Nx, Ny)

    ds.x = ib.dxi';
    ds.y = ib.deta';
    
    beta = hypot(ds.x, ds.y);
    
    %% u-velocity
%     [~, xi_x] = min(abs(ib.xi - grid.Xu), [], 2)
%     [~, eta_x] = min(abs(ib.eta - grid.Yu), [], 2);
    
    d.x = zeros(size(grid.dXv));

    for i = 1:length(ib.xi)
        
        r.x = ib.xi(i) - grid.Xu;
        r.y = ib.eta(i) - grid.Yu;
        
        d.x(i,:) = delta(r.x, grid.dXv);
        d.y(i,:) = delta(r.y, grid.dY');
        
    end

    dr.x = repmat(grid.dXv, [Ny-1,1]);
    dr.y = repmat(grid.dY, [Nx,1]);

    H.u = beta.*delta(r.x, dr.x).*delta(r.y, dr.y);

    %% v-velocity
    r.x = ib.xi' - reshape(grid.xv, [], 1);
    r.y = ib.eta' - reshape(grid.yv, [], 1);
    
    dr.x = repmat(grid.dXv, [1,Ny])';
    dr.y = repmat(grid.dYu', [1,Nx-1])';
    
    H.v = beta.*delta(r.x, dr.x).*delta(r.y, dr.y);

end

