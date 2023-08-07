function E = Ehat(grid, ib, Nx, Ny)
%% u-velocity

Eux = interpolation1D(ib.xi, grid.Xu, grid.dX, 1);
Euy = interpolation1D(ib.eta, grid.Yu, grid.dYp, 1);

Euyy = sparse(kron(ones(size(grid.Xu)), Euy));
Euxx = sparse(kron(Eux, ones(size(grid.Yu))));

E.u = Euyy.*Euxx;
    
%% v-velocity
    
Evx = interpolation1D(ib.xi, grid.Xv, grid.dXp, 1);
Evy = interpolation1D(ib.eta, grid.Yv, grid.dY, 1);

Evyy = sparse(kron(ones(size(grid.Xv)), Evy));
Evxx = sparse(kron(Evx, ones(size(grid.Yv))));

E.v = Evyy.*Evxx;

%% Matrix assembly

% E = blkdiag([E.u E.v]);

end

