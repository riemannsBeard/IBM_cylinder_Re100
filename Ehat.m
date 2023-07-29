function E = Ehat(grid, ib, Nx, Ny)


    %% u-velocity
%     r.x = reshape(grid.X - ib.xi, [], 1);
%     r.y = reshape(grid.Y - ib.eta, [], 1);

for i = 1:length(ib.xi)
    
    %% u-velocity
    
    
%     r.x = reshape(grid.xv, 1, []) - ib.xi(i);
%     r.y = reshape(grid.yv, 1, []) - ib.eta(i);

    r.x = grid.xu - ib.xi(i);
    r.y = grid.yu - ib.eta(i);

    
    dr.x = repmat(grid.dXv, Nx, 1);
    dr.y = repmat(grid.dY', Ny-1, 1)';

%     dr.x = repmat(grid.dX, [1,numel(ib.xi)]);
%     dr.y = repmat(grid.dY, [1,numel(ib.xi)]);
    
    alpha.u = dr.x.*dr.y;

    E.u(i,:) = reshape(alpha.u.*delta(r.x, dr.x).*delta(r.y, dr.y),[], 1);
    
    %% v-velocity

    r.x = grid.xv - ib.xi(i);
    r.y = grid.yv - ib.eta(i);
    
    dr.x = repmat(grid.dX', Nx-1, 1);
    dr.y = repmat(grid.dYu, Ny, 1)';

%     dr.x = repmat(grid.dX, [1,numel(ib.xi)]);
%     dr.y = repmat(grid.dY, [1,numel(ib.xi)]);
    
    alpha.v = dr.x.*dr.y;

    E.v(i,:) = reshape(alpha.v.*delta(r.x, dr.x).*delta(r.y, dr.y),[], 1);
    

end
    
    %% v-velocity
%     r.x = reshape(grid.xv, 1, []) - ib.xi;
%     r.y = reshape(grid.yv, 1, []) - ib.eta;
%     
%     dr.x = repmat(grid.dXv', [Ny,1])';
%     dr.y = repmat(grid.dY, [Nx-1,1])';
%     
%     alpha.v = dr.x.*dr.y;
% 
%     E.v = alpha.v.*delta(r.x, dr.x).*delta(r.y, dr.y);
    


end

