function E = Ehat(grid, ib, Nx, Ny)

% % u-velocity
% 
% Eux = interpolation1D(ib.xi, grid.Xu, grid.dX, 1);
% Euy = interpolation1D(ib.eta, grid.Yu, grid.dYp, 1);
% 
% Euyy = kron(ones(size(grid.Xu)), sparse(Euy));
% Euxx = kron(sparse(Eux), ones(size(grid.Yu)));
% E.u = Euyy.*Euxx;
% 
% % v-velocity
% 
% Evx = interpolation1D(ib.xi, grid.Xv, grid.dXp, 1);
% Evy = interpolation1D(ib.eta, grid.Yv, grid.dY, 1);
% 
% Evyy = kron(ones(size(grid.Xv)), sparse(Evy));
% Evxx = kron(sparse(Evx), ones(size(grid.Yv)));
% E.v = Evyy.*Evxx;

    %% u-velocity
    
    dr.x = repmat(grid.dX', [1,Ny-1]);
    dr.y = repmat(grid.dYu, [1,Nx]);
    
    r.x = reshape(grid.xu, 1, []) - ib.xi;
    r.y = reshape(grid.yu, 1, []) - ib.eta;
       
    E.u = dr.x.*dr.y.*delta(r.x, dr.x).*delta(r.y, dr.y);

    %% v-velocity

    dr.x = repmat(grid.dXv', [Ny,1])';
    dr.y = repmat(grid.dY, [Nx-1,1])';    
    
    r.x = reshape(grid.xv, 1, []) - ib.xi;
    r.y = reshape(grid.yv, 1, []) - ib.eta;
    
    E.v = dr.x.*dr.y.*delta(r.x, dr.x).*delta(r.y, dr.y);

end

