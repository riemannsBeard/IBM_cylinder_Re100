function [H, beta] = Hhat(grid, ib, Nx, Ny)
%% u-velocity

Hux = interpolation1D(ib.xi, grid.Xu, grid.dX, 0);
Huy = interpolation1D(ib.eta, grid.Yu, grid.dYp, 0);

Huyy = sparse(kron(ones(size(grid.Xu)), Huy));
Huxx = sparse(kron(Hux, ones(size(grid.Yu))));

H.u = Huxx.*Huyy;
    
%% v-velocity
    
Hvx = interpolation1D(ib.xi, grid.Xv, grid.dXp, 0);
Hvy = interpolation1D(ib.eta, grid.Yv, grid.dY, 0);

Hvyy = sparse(kron(ones(size(grid.Xv)), Hvy));
Hvxx = sparse(kron(Hvx, ones(size(grid.Yv))));

H.v = Hvxx.*Hvyy;

%% Scaling factor beta

beta = sparse(blkdiag(diag(ib.ds, 0), diag(ib.ds, 0)));

%% Matrix assembly

% H = blkdiag([H.u H.v]);

end

