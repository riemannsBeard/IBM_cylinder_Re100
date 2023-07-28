function [H] = Hhat(grid, ib, Nx, Ny)
    
%     % u-velocity
% 
%     Hux = interpolation1D(ib.xi, grid.Xu, grid.dX, 0);
%     Huy = interpolation1D(ib.eta, grid.Yu, grid.dYp, 0);
%     
%     Huyy = kron(ones(size(grid.Xu)), sparse(Huy));
%     Huxx = kron(sparse(Hux), ones(size(grid.Yu)));
%     H.u = Huyy.*Huxx;
%     
%     % v-velocity
%         
%     Hvx = interpolation1D(ib.xi, grid.Xv, grid.dXp, 0);
%     Hvy = interpolation1D(ib.eta, grid.Yv, grid.dY, 0);
%     
%     Hvyy = kron(ones(size(grid.Xv)), sparse(Hvy));
%     Hvxx = kron(sparse(Hvx), ones(size(grid.Yv)));
%     H.v = Hvyy.*Hvxx;
%     
%     diag_ds = spdiags(ib.ds(:), 0, numel(ib.ds), numel(ib.ds));
%     
%     beta = blkdiag(diag_ds, diag_ds);


    ds.x = ib.dxi';
    ds.y = ib.deta';
    
    beta = hypot(ds.x, ds.y);
    
    %% u-velocity
    r.x = ib.xi' - reshape(grid.xu, [], 1);
    r.y = ib.eta' - reshape(grid.yu, [], 1);
    
    dr.x = repmat(grid.dX, [Ny-1,1]);
    dr.y = repmat(grid.dYu', [Nx,1]);

    H.u = beta.*delta(r.x, dr.x).*delta(r.y, dr.y);

    %% v-velocity
    r.x = ib.xi' - reshape(grid.xv, [], 1);
    r.y = ib.eta' - reshape(grid.yv, [], 1);
    
    dr.x = repmat(grid.dXv, [1,Ny])';
    dr.y = repmat(grid.dY', [1,Nx-1])';
    
    H.v = beta.*delta(r.x, dr.x).*delta(r.y, dr.y);    
    
end

