clc
clear
close all

%% Defaults
set(0, 'DefaultLineLineWidth', 2);
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultLegendFontSize', 14);
set(0, 'DefaultLegendFontSizeMode', 'manual');
set(0, 'DefaultAxesTickLabelInterpreter', 'latex');
set(0, 'DefaultLegendInterpreter', 'latex');
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'DefaultAxesFontsize', 14);
set(0, 'DefaultAxesLineWidth', 1.0);

%% Datos

Re = 100;
Nx = 1024; % Celdillas en x
Ny = 1024; % Celdillas en y
Lx = 26;
Ly = 26;
tSampling = 0.5;
CFL = 0.5;

%% Staggered Grid Generation

alpha.x = 0;
alpha.y = 0;

[grid, u, v, p] = gridGeneration(Lx, Ly, Nx, Ny, alpha);
%dt = CFL*min(grid.cellMin^2*Re, grid.cellMin);
dt = 1e-3;

itSampling = floor(tSampling/dt);

%% Immersed Boundary Grid Generation

Nk = 128;
R = 0.5;
ib = IBMcylinder(R, Nk);
ib.xi = ib.xi + 6.5;
ib.eta = ib.eta + 0.5*Ly;

%% Boundary conditions

bc.uS = ones(1,Nx-1);
bc.uN = ones(1,Nx-1);
bc.uE = zeros(Ny,1);
bc.uW = ones(Ny,1);

bc.vS = zeros(1,Nx);
bc.vN = zeros(1,Nx);
bc.vE = zeros(Ny-1,1);
bc.vW = zeros(Ny-1,1);

%% Immersed Boundary Conditions

ib.Utheta = 0*pi/8;

ib.u = ib.Utheta*sin(ib.theta);
ib.v = -ib.Utheta*cos(ib.theta);

%% Grid detail
% lx = 2;
% l = grid.x < lx;
% fig = figure;
% plot(grid.x(l), grid.y(l), 'k.'), pbaspect([lx Ly 1])
% set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
% xlabel('$l_x$', 'interpreter', 'latex', 'fontsize', 16)
% ylabel('$L_y$', 'interpreter', 'latex', 'fontsize', 16)
% printFigure, print(fig, 'pipeFlow_gridDetail', '-dpdf', '-r0');


%% Operators

[D, G, R, M, bc2_] = DGRM(grid, Nx, Ny);

Lhat = laplacianHat(Nx, Ny, grid);
L.L = M.hat*Lhat.L/R.R;

M.M = M.hat/R.R;
M.inv = inv(M.M);

Ahat = sparse(speye(size(Lhat.L))/dt - 0.5*Lhat.L/Re);
A = M.hat*Ahat/R.R;

BN = dt*speye(size(M.M))/M.M + (0.5/Re)*dt*dt*(M.inv*L.L)*M.inv +...
    ((0.5/Re)^2)*(dt^3)*((M.inv*L.L)^2)*M.inv;

% BC's due to Laplacian
bc1hat.u = Lhat.ux0*bc.uW + Lhat.uy1*bc.uN' + ...
    Lhat.uy0*bc.uS';
bc1hat.v = Lhat.vx0*bc.vW + Lhat.vx1*bc.vE + Lhat.vy1*bc.vN' + ...
    Lhat.vy0*bc.vS';

bc1 = M.hat*[bc1hat.u; bc1hat.v]/Re;

% BC's due to Divergence
bc2 = D.uW*(bc.uW.*grid.dY) + ...
    D.vS*(bc.vS'.*grid.dX) + D.vN*(bc.vN'.*grid.dX);

r2 = [-bc2; ib.u; ib.v];

% LU decomposition for performance
% A = decomposition(A);
% R.R = decomposition(R.R);

%% IBM stuff

% Regularization
[Hhat_, beta] = Hhat(grid, ib, Nx, Ny);
Hhat_ = blkdiag(Hhat_.u, Hhat_.v);
H = sparse(M.M*Hhat_);

% Interpolation
[Ehat_, alpha] = Ehat(grid, ib, Nx, Ny);
Ehat_ = blkdiag(Ehat_.u, Ehat_.v);
E = sparse(Ehat_/R.R);

% LU decomposition for performance
EH = sparse(E*H);
% EH = decomposition(EH);

%% Left-Hand Side term

Q = [G.G, E'];

LHS = sparse(Q'*BN*Q);
% LHS = decomposition(LHS);

%% Simulation

u = reshape(u, [], 1);
v = reshape(v, [], 1);

uOld = u;
vOld = v;

t0 = 0;
tf = 150;
t = dt:dt:tf;

% Preallocation for effieciency
epsU = t*0;
epsV = t*0;
F = t*0;
f.x = t*0;
f.y = t*0;

% Matrix product to increase performance
EHEE = (EH)\E*E';
BNQ = BN*Q;

%% Temporal evolution
tic
for k = 1:length(t)
            
    % Advective terms
    [NhatOld, ~, ~] = convectionHat(grid, uOld, vOld, Nx, Ny, bc);
    [Nhat, ua, va] = convectionHat(grid, u, v, Nx, Ny, bc);
    
    rnHat = explicitTerms(Lhat, Re, dt, Nhat, NhatOld, u, v);  
    rn = M.hat*rnHat;
        
    %% 1. Solve for intermediate velocity
       
    r1 = rn + bc1;
    
    % Flux calculation    
    q = A\r1;
    qu = q(1:Ny*(Nx-1));
    qv = q(Ny*(Nx-1)+1:end);

    %% 2. Solve the Poisson Equation
    
    RHS = Q'*q - r2;
    lambda = LHS\RHS;
    
    % Pressure calculation
    phi = lambda(1:end-2*Nk);
    
    % Forces
    fTilda.x = lambda(end-2*Nk+1:end-Nk);
    fTilda.y = lambda(end-Nk+1:end);
    fTilda.f = [fTilda.x; fTilda.y];
    
%     f.f = -(EH)\E*E'*fTilda.f;
    f.f = -EHEE*fTilda.f;    
    
    %% 3. Projection step
    
%     q = q - BN*Q*lambda;
    q = q - BNQ*lambda;

    vel = R.R\q;

    % Residuals
    epsU(k) = max(abs(u - vel(1:Ny*(Nx-1))));
    epsV(k) = max(abs(v - vel(Ny*(Nx-1)+1:end)));

    % Separation of velocity components
    u = vel(1:Ny*(Nx-1));
    v = vel(Ny*(Nx-1)+1:end);
    phi = reshape(phi, Ny, Nx);

    % Forces storage
    f.x(k) = sum(f.f(1:Nk));
    f.y(k) = sum(f.f(Nk+1:end));
    F(k) = hypot(f.x(k), f.y(k));    
    
    % Information
    disp(['t = ' num2str(t(k))]); % '. Elapsed time: ' num2str(toc) 's.']);
    disp(['Residuals u = ' num2str(epsU(k))]);
    disp(['Residuals v = ' num2str(epsV(k)) newline]);

    % On-the-fly plots (comment for speeding up the code)
%     if (mod(k, itSampling) == 0)
%         figure(1),
%         
%         ax = subplot(311);
%         pcolor(grid.x, grid.y, ua), shading interp
%         cmap = jet;
%         colorbar, colormap(ax, cmap)
%         set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
%         title('$u$', 'interpreter', 'latex', 'fontsize', 16)
%         pbaspect([Lx Ly 1])
%         
%         ax = subplot(312);
%         pcolor(grid.x, grid.y, va), shading interp
%         cmap = jet;
%         colorbar, colormap(ax, cmap)
%         set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
%         title('$v$', 'interpreter', 'latex', 'fontsize', 16)
%         pbaspect([Lx Ly 1])
%         
%         ax = subplot(313);
%         contourf(grid.xp, grid.yp, phi), %shading interp,
%         colorbar
% %         cmap = jet;
% %         colorbar, colormap(ax, cmap)
%         set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
%         title('$\phi$', 'interpreter', 'latex', 'fontsize', 16)
%         pbaspect([Lx Ly 1])
%         
%         drawnow
%         
%         figure(2),
%         loglog(1:k, epsU(1:k), 1:k, epsV(1:k))
%         set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
%         h = legend('$u$', '$v$');
%         set(h, 'interpreter', 'latex', 'fontsize', 16)
%         xlabel('$N$', 'interpreter', 'latex', 'fontsize', 16)
%         ylabel('$\xi$', 'interpreter', 'latex', 'fontsize', 16)
%         title('Residuals')
%         drawnow
%         
%     end
    
end
toc

%% Plots

% Contours
fig = figure(1);
subplot(211)
contourf(grid.x - 6.5, grid.y - Ly/2, hypot(ua,va), 32,...
    'LineStyle', 'none'),
hold on
shading interp,
xlim([-5 18]), ylim([-4. 4.])
plot(ib.xi - 6.5, ib.eta - Ly/2, 'w-', 'linewidth', 1.25)
box on
colormap('jet');
c = colorbar;
c.TickLabelInterpreter = 'latex';
title('$\mathbf{u}$', 'interpreter', 'latex', 'fontsize', 16)
% pbaspect([Lx Ly 1])
saveas(fig, 'cylinder_Re100', 'jpeg');


% vorticity = (diff(v'))./(diff(grid.xv')) - diff(u)./diff(grid.yu);

% % Contours
% fig = figure(2);
% subplot(211)
% contourf(grid.x - 6.5, grid.y - Ly/2, vorticity, 32,...
%     'LineStyle', 'none'),
% hold on
% shading interp,
% xlim([-5 18]), ylim([-4. 4.])
% plot(ib.xi - 6.5, ib.eta - Ly/2, 'w-', 'linewidth', 1.25)
% box on
% colormap('jet');
% c = colorbar;
% c.TickLabelInterpreter = 'latex';
% title('$\mathbf{u}$', 'interpreter', 'latex', 'fontsize', 16)
% % pbaspect([Lx Ly 1])
% saveas(fig, 'cylinder_Re100', 'jpeg');

% Forces
figure,
plot(t, 2*f.y, t, -2*f.x)
xlabel('$\tau$')
legend('$C_L$', '$C_D$', 'interpreter', 'latex')

% Residuals
fig = figure(5);
loglog(1:k, epsU(1:k), 1:k, epsV(1:k))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
h = legend('$u$', '$v$');
set(h, 'interpreter', 'latex', 'fontsize', 16)
xlabel('$N$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\xi$', 'interpreter', 'latex', 'fontsize', 16)
title('Residuals')
printFigure, print(fig, 'cylinder_Re100_res', '-dpdf', '-r0');



