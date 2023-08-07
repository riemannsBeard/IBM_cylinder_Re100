function [ grid, u, v, p ] = gridGeneration(Lx, Ly, Nx, Ny, hmin, x0)

%     grid.X = linspace(0, Lx, Nx+1);
%     grid.Y = linspace(0, Ly, Ny+1);

%     r = (30/(x0 - hmin))^(1/(131 - 1));
%     x = (x0 + hmin)*r.^(0:(131 - 1));
%     grid.X = [-flip(x), -x0:hmin:x0, x] + 30; grid.X(1) = 0;
%     grid.Y = grid.X;


% Reference: "A new compact difference scheme for second derivative in non-
% uniform grid expressed in self-adjoint form" by TK Sengupta, S Bhaumik 
% and Shameem Usman, J. Comp. Phys. 230 (2011) 1822–1848
    N = 110;
    beta = 2.25;
 
    x = (0.5*Lx - x0)*(1 - tanh(beta*(1 - (1:N)/N))./tanh(beta));
    xr = x + 0.5*Lx + x0; xr(end) = Lx;
    xl = flip(0.5*Lx - x) - x0; xl(1) = 0;
    grid.X = [xl, (0.5*Lx - x0):hmin:(x0 + 0.5*Lx), xr];
    grid.X = grid.X - 0.5*Lx;
    grid.Y = grid.X;

    grid.dX = diff(grid.X)';
    grid.dY = diff(grid.Y)';

    grid.cellMin = min([grid.dX; grid.dY]);

    grid.Xp = 0.5*(grid.X(2:end) + grid.X(1:end-1));
    grid.Yp = 0.5*(grid.Y(2:end) + grid.Y(1:end-1));
    
    grid.dXp = diff(grid.Xp)';
    grid.dYp = diff(grid.Yp)';
    
    p = zeros(Ny, Nx);
    u = ones(Ny, Nx-1);
    v = zeros(Ny-1, Nx);

    grid.Xu = grid.X(2:end-1);
    grid.Yu = grid.Yp;
    
    grid.dXu = diff(grid.Xu);
    grid.dYu = diff(grid.Yu);      

    grid.Xv = grid.Xp;
    grid.Yv = grid.Y(2:end-1);
    
    grid.dXv = diff(grid.Xv);
    grid.dYv = diff(grid.Yv);    

    [grid.x, grid.y] = meshgrid(grid.X, grid.Y);
    [grid.xp, grid.yp] = meshgrid(grid.Xp, grid.Yp);
    [grid.xu, grid.yu] = meshgrid(grid.Xu, grid.Yu);
    [grid.xv, grid.yv] = meshgrid(grid.Xv, grid.Yv);
    
    grid.dx = diff(grid.x')';
    grid.dy = diff(grid.y);

    grid.dxp = diff(grid.xp')';
    grid.dyp = diff(grid.yp);
    
end

