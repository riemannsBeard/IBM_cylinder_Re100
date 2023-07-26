function [ D, G, R, M ] = DGRM( grid, Nx, Ny )
% Funcion que devuelve los operadores matriciales G, D, R y M
    %% OPERADOR GRADIENTE
    % Se aplica a la presion
        
    % Operador derivada para una fila y una columna
    g.x = spdiags(ones(Nx-1, 2).*[-1, 1], [0, 1], Nx-1, Nx);
    g.y = spdiags(ones(Ny-1, 2).*[-1, 1], [0, 1], Ny-1, Ny);
           
    %% OPERADOR DIVERGENCIA
%     % Se aplica a las velocidades
%     d.x = -g.x';
%     d.y = -g.y';
%     
%     % Account for Neumann BC
%     d.x(end,end) = 0;
% 
    G.x = kron(speye(Ny), g.x);
    G.y = kron(g.y, speye(Nx));
% 
%     D.D = [D.x D.y];
    
    G.G = [G.x; G.y];
    
%     bcW = zeros(Nx,1); bcW(1) = 1;
%     bcE = zeros(Nx,1); bcE(end) = -1;
%     bcN = zeros(Ny,1); bcN(1) = 1;
%     bcS = zeros(Ny,1); bcS(end) = -1;
%     
%     bc2.uW = kron(bcW, grid.dY);
%     bc2.uE = kron(grid.dY, bcE);
%     bc2.vN = kron(bcN, grid.dX);
%     bc2.vS = kron(bcS, grid.dY);
    
    % Boundary terms
    D.uW = kron(speye(Ny), spdiags(-1, 0, Nx, 1));
    D.uE = kron(speye(Ny), spdiags(1, -Nx+1, Nx, 1));
    D.vS = kron(spdiags(-1, 0, Ny, 1), speye(Nx));
    D.vN = kron(spdiags(1, -Ny+1, Ny, 1), speye(Nx));
        
    %% MATRIZ DE FLUJO
    
    dyj = sparse(diag(grid.dY, 0));
    dxi = sparse(diag(grid.dX, 0));
    
    R.u = kron(dyj, speye(Nx-1));
    R.v = kron(speye(Ny-1), dxi);

    R = blkdiag(R.u, R.v);
           
    %% MATRIZ DE MASA
    
    ix = [0.75; ones(Nx-2,1); 0.75];
    iy = [0.75; ones(Ny-2,1); 0.75];

    Ix = sparse(diag(ix, 0));
    Iy = sparse(diag(iy, 0));
    
    Dxp = sparse(diag(grid.dXp, 0));
    Dyp = sparse(diag(grid.dYp, 0));

    % Mhat
    Mhat.u = kron(Iy, Dxp);
    Mhat.v = kron(Dyp, Ix);

    M.hat = blkdiag(Mhat.u, Mhat.v);
    
    % M
    M.M = R\M.hat;
    
end

