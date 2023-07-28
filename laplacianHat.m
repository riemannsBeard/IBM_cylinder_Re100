function [ L ] = laplacianHat( Nx, Ny, grid )
%     % Vectores unitarios auxiliares
%     ex = ones(Nx+2, 1);
%     ey = ones(Ny+2, 1);
    
    %% Laplaciano de u
        
    % Lux
    gux = spdiags([-1./grid.dX 1./grid.dX], [0 1], Nx, Nx+1);
    guxx = spdiags([-1./grid.dXp 1./grid.dXp], [0 1], Nx-1, Nx)*gux;
    L.ux = kron(guxx*spdiags(ones(Nx-1, 1), -1, Nx+1, Nx-1), speye(Ny));
    L.ux0 = kron(guxx*spdiags(1, 0, Nx+1, 1), speye(Ny));
    L.ux1 = kron(guxx*spdiags(1, -Nx, Nx+1, 1), speye(Ny));
    
       
    % Luy
    dYp_ = [0.5*grid.dY(1); grid.dYp; 0.5*grid.dY(end)];
    guy = spdiags([-1./dYp_ 1./dYp_], [0 1], Ny+1, Ny+2);
    dY_ = [0.75*grid.dY(1); grid.dY(2:end-1); 0.75*grid.dY(end)];
    guyy = spdiags([-1./dY_ 1./dY_], [0 1], Ny, Ny+1)*guy;
    
    L.uy = kron(speye(Nx-1), guyy*spdiags(ones(Ny, 1), -1, Ny+2, Ny));
    L.uy0 = kron(speye(Nx-1), guyy*spdiags(1, 0, Ny+2, 1));
    L.uy1 = kron(speye(Nx-1), guyy*spdiags(1, -(Ny+1), Ny+2, 1));
    
    % Operador Laplaciano de u
    L.u = L.ux + L.uy;
    
    %% Laplaciano de v
    
    % Lvx
    dXp_ = [0.5*grid.dX(1); grid.dXp; 0.5*grid.dX(end)];
    gvx = spdiags([-1./dXp_ 1./dXp_], [0 1], Nx+1, Nx+2);
    dX_ = [0.75*grid.dX(1); grid.dX(2:end-1); 0.75*grid.dX(end)];
    gvxx = spdiags([-1./dX_ 1./dX_], [0 1], Nx, Nx+1)*gvx;
        
    L.vx = kron(gvxx*spdiags(ones(Nx, 1), -1, Nx+2, Nx), speye(Ny-1));
    L.vx0 = kron(gvxx*spdiags(1, 0, Nx+2, 1), speye(Ny-1));
    L.vx1 = kron(gvxx*spdiags(1, -(Nx+1), Nx+2, 1), speye(Ny-1));    
    
    % Lvy
    gvy = spdiags([-1./grid.dY 1./grid.dY], [0 1], Ny, Ny+1);
    gvyy = spdiags([-1./grid.dYp 1./grid.dYp], [0 1], Ny-1, Ny)*gvy;
    
    L.vy = kron(speye(Nx), gvyy*spdiags(ones(Ny-1, 1), -1, Ny+1, Ny-1));
    L.vy0 = kron(speye(Nx), gvyy*spdiags(1, 0, Ny+1, 1));
    L.vy1 = kron(speye(Nx), gvyy*spdiags(1, -Ny, Ny+1, 1));
        
    % Operador Laplaciano de u
    L.v = L.vx + L.vy;

%%    
    L.L = blkdiag(L.u, L.v);
        
end

