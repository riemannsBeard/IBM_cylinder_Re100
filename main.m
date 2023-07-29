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

%% Load/build matrix stuff
if ~exist('./matrixStuff.mat', 'file')
    %% Datos
    Re = 200;
    Nx = 300; % Celdillas en x
    Ny = 300; % Celdillas en y
    hmin = 0.02;
    x0 = 0.8;
    Lx = 60;
    Ly = 60;
    tSave = 1;
    CFL = 0.4;
    
    %% Staggered Grid Generation
    [grid, u, v, p] = gridGeneration(Lx, Ly, Nx, Ny, hmin, x0);
    %dt = CFL*min(grid.cellMin^2*Re, grid.cellMin);
    dt = 5e-3;

%     dt1 = CFL*min(grid.dX);
%     dt2 = CFL*min(grid.dY);
%     dt = min(dt1,dt2);
    
    %itSampling = floor(tSampling/dt);
    
    %% Immersed Boundary Grid Generation
    r = 0.5;
    Nk = round(2*pi*r/grid.cellMin);

    ib = IBMcylinder(r, Nk);
    ib.xi = ib.xi + 0.5*Lx;
    ib.eta = ib.eta + 0.5*Ly;
    ib.ds = 2*pi*r./Nk.*ones(Nk, 1);
    
    %% Boundary conditions
    bc.uS = ones(1,Nx-1);
    bc.uN = ones(1,Nx-1);
    bc.uE = ones(Ny,1);
    bc.uW = ones(Ny,1);
    
    bc.vS = zeros(1,Nx);
    bc.vN = zeros(1,Nx);
    bc.vE = zeros(Ny-1,1);
    bc.vW = zeros(Ny-1,1);
    
    %% Immersed Boundary Conditions
    ib.Utheta = 0*pi/8;
    
    ib.u = 0*ib.Utheta*sin(ib.theta);
    ib.v = 0*(-ib.Utheta*cos(ib.theta));
    
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
    [D, G, R, M] = DGRM(grid, Nx, Ny);
    
    Lhat = laplacianHat(Nx, Ny, grid);
    L.L = M.hat*Lhat.L/R;
    
    M.M = M.hat/R;
    M.inv = inv(M.M);
    
    Ahat = sparse(speye(size(Lhat.L))/dt - 0.5*Lhat.L/Re);
    A = M.hat*Ahat/R;
    clear Ahat
    
    BN = dt*speye(size(M.M))/M.M + (0.5/Re)*dt*dt*(M.inv*L.L)*M.inv +...
        ((0.5/Re)^2)*(dt^3)*((M.inv*L.L)^2)*M.inv;
    clear L
    
%    % BC's due to Laplacian
%     bc1hat.u = Lhat.ux0*bc.uW + Lhat.uy1*bc.uN' + ...
%         Lhat.uy0*bc.uS';
%     bc1hat.v = Lhat.vx0*bc.vW + Lhat.vx1*bc.vE + Lhat.vy1*bc.vN' + ...
%         Lhat.vy0*bc.vS';
%     
%     bc1 = M.hat*[bc1hat.u; bc1hat.v]/Re;
%     clear bc1hat
    
%     % BC's due to Divergence
%     bc2 = -(D.uW*(bc.uW.*grid.dY) + D.uE*(r.*grid.dY) + ...
%         D.vS*(bc.vS'.*grid.dX) + D.vN*(bc.vN'.*grid.dX));
%     
%     r2 = [-bc2; ib.u; ib.v];
%     clear D
    
    %% IBM stuff
    % Regularization
    [Hhat_, beta] = Hhat(grid, ib, Nx, Ny);
    Hhat_ = blkdiag(Hhat_.u, Hhat_.v);
    H = sparse(M.M*Hhat_);
    Mhat = M.hat;
    clear Hhat_
    
    % Interpolation
    Ehat_ = Ehat(grid, ib, Nx, Ny);
    Ehat_ = blkdiag(Ehat_.u, Ehat_.v);
    E = sparse(Ehat_/R);
    clear Ehat_
    
    H = E';
    % Matrix product to increase performance
    EH = sparse(E*H);
    clear H
    
    % Matrix product to increase performance
    EHEE = (EH)\E*E';
    clear EH
    
    %% Left-Hand Side term
    Q = [G.G, E'];
    clear G E
    
    % Matrix product to increase performance
    BNQ = BN*Q;
    LHS = sparse(Q'*BN*Q);
    clear BN
    
    % Save matrix stuff 
    save('./matrixStuff.mat')
    
else
    % Load matrix stuff
    load('./matrixStuff.mat')
    
end


%% Simulation
u = reshape(u, [], 1);
v = reshape(v, [], 1);

uOld = u;
vOld = v;

t0 = 0;
tf = 5;
t = t0:dt:tf;

% Preallocation for effieciency
epsR = NaN(size(t));
epsU = NaN(size(t));
epsV = NaN(size(t));
F = NaN(size(t));
f.x = NaN(size(t));
f.y = NaN(size(t));

%% Create storage and forces dir
[~, ~, ~] = mkdir('stored');
[~, ~, ~] = mkdir('forces');

% Initialize forces evolution file
forcesID = fopen('./forces/forces0', 'w');
fprintf(forcesID, 'time \t fx \t fy\n');

q = R*[u(:); v(:)];
qast = q;

% Advective terms
[NhatOld, ~, ~, ~, ~] = advectionHat(grid, uOld, vOld, Nx, Ny, bc);
% [Nhat, ua, va, ue, ve] = advectionHat(grid, u, v, Nx, Ny, bc);
Nhat = NhatOld;

%% LU decomposition (not allowed to be saved)
LHS = decomposition(LHS);
A = decomposition(A);
R = decomposition(R);

%% Temporal evolution

tic
for k = 1:length(t)
            
    u = reshape(u, [], 1);
        
    bc1hat.u = Lhat.ux0*bc.uW + Lhat.ux1*bc.uE + Lhat.uy1*bc.uN' + ...
        Lhat.uy0*bc.uS';
    bc1hat.v = Lhat.vx0*bc.vW + Lhat.vx1*bc.vE + Lhat.vy1*bc.vN' + ...
        Lhat.vy0*bc.vS';
    
    bc1 = M.hat*[bc1hat.u; bc1hat.v]/Re;
    
    rnHat = explicitTerms(Lhat, Re, dt, Nhat, NhatOld, u, v);  
    rn = Mhat*rnHat;
        
    %% 1. Solve for intermediate velocity
    r1 = rn + bc1;
    
    % Flux calculation    
    qast = A\r1;
    qu = q(1:Ny*(Nx-1));
    qv = q(Ny*(Nx-1)+1:end);

    %% 2. Solve the Poisson Equation
    % BC's due to Divergence
    bc2 = -(D.uW*(bc.uW.*grid.dY) + D.uE*(bc.uE.*grid.dY) + ...
        D.vS*(bc.vS'.*grid.dX) + D.vN*(bc.vN'.*grid.dX));
    
    r2 = [-bc2; ib.u; ib.v];    
    
    RHS = Q'*qast - r2;
    lambda = LHS\RHS;
        
    % Pressure calculation
    phi = lambda(1:end-2*Nk);
    
    % Forces
    fTilda.x = lambda(end-2*Nk+1:end-Nk);
    fTilda.y = lambda(end-Nk+1:end);
    fTilda.f = [fTilda.x; fTilda.y];
    
    f.f = -EHEE*fTilda.f;
    
    %% 3. Projection step
    % Flux and velocity calculation
    qp1 = qast - BNQ*lambda;
    
    % Residuals calculation
    epsR(k) = norm(qp1 - q)/(dt*norm(qp1));

    % Forces storage
    f.x(k) = sum(f.f(1:Nk));
    f.y(k) = sum(f.f(Nk+1:end));
    F(k) = hypot(f.x(k), f.y(k));
    
    % Update flux
    q = qp1;
    
    % Recover velocity
    vel = R\q;
    
%     % Recover pressure
%     phi = reshape(phi, Ny, Nx);    
    
    % Residuals
    epsU(k) = norm(u - vel(1:Ny*(Nx-1)));
    epsV(k) = norm(v - vel(Ny*(Nx-1)+1:end));    
    
    % Separation of velocity components
    u = vel(1:Ny*(Nx-1));
    v = vel(Ny*(Nx-1)+1:end);
    
    u = reshape(u, Ny, Nx-1);
    
    bc.uE = bc.uE - dt./grid.dX.*(bc.uE - u(:,end));

    % Advective terms
    NhatOld = Nhat; 
    [Nhat, ua, va, ue, ve] = advectionHat(grid, u, v, Nx, Ny, bc);

%     u = reshape(u, Ny, Nx-1);
%     bc.uE(:,end) = u(:,end) - dt*(u(:,end) - u(:,end-1))./grid
% 
%     UE = u[:,-1] - Uinf*(u[:,-1]-u[:,-2])*dt/np.diff(Xu)[0,-1]

    % Forces writing
    fprintf(forcesID, '\n%6.6f \t %6.6f \t %6.6f', [t(k) -f.x(k) -f.y(k)]);
    
    % Information
    disp(['t = ' num2str(t(k))]); % '. Elapsed time: ' num2str(toc) 's.']);
    disp(['Residuals q = ' num2str(epsR(k))]);    
    disp(['Residuals u = ' num2str(epsU(k))]);
    disp(['Residuals v = ' num2str(epsV(k)) newline]);
    
    if t(k) > (tSave - dt) && t(k) <= (tSave + dt)
                
        vorticity = computeVorticity(ua, va, grid);
        
        save(['./stored/cylinder_Re100_t' num2str(t(k)) '.mat'],...
            'ua', 'va', 'vorticity', 'grid');
        
    end

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
% %         cmap = jet;vim
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

fclose(forcesID);

%% Plots

% Contours
fig = figure(1);
subplot(211)
contourf(grid.x, grid.y, hypot(ua,va), 32,...
    'LineStyle', 'none'),hold on
% plot(grid.x, grid.y, 'k.')
xlabel('$x$'), ylabel('$y$')
hold on
shading interp,
xlim(0.5*Lx + [-7 14]), ylim(0.5*Ly + [-5 5])
plot(ib.xi, ib.eta, 'w-', 'linewidth', 1.25)
box on
colormap('jet');
c = colorbar;
c.TickLabelInterpreter = 'latex';
c.Label.Interpreter = 'latex';
c.Label.String = '$\mathbf{u}$';
c.Label.FontSize = 16;
% saveas(fig, 'cylinder_Re100', 'jpeg');

%%
vorticity = computeVorticity(ua, va, grid);
subplot(212)
contourf(grid.x, grid.y, vorticity, 16,...
    'LineStyle', 'none'),
xlabel('$x$'), ylabel('$y$')
hold on
xlim(0.5*Lx + [-7 14]), ylim(0.5*Ly + [-5 5])
plot(ib.xi, ib.eta, 'w-', 'linewidth', 1.25)
box on
colormap('jet');
c = colorbar;
caxis([-8 8])
c.TickLabelInterpreter = 'latex';
c.Label.Interpreter = 'latex';
c.Label.String = '$\omega_z$';
c.Label.FontSize = 16;
saveas(fig, 'cylinder_Re100', 'jpeg');

%% Forces
figure,
plot(t, -2*f.y, t, -2*f.x)
ylim([-0.5 5])
xlabel('$\tau$')
legend('$C_L$', '$C_D$', 'interpreter', 'latex')

%% Residuals
fig = figure(5);
loglog(1:k, epsU(1:k), 1:k, epsV(1:k))
set(gca, 'TickLabelInterpreter','latex', 'fontsize', 12)
h = legend('$u$', '$v$');
set(h, 'interpreter', 'latex', 'fontsize', 16)
xlabel('$N$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\xi$', 'interpreter', 'latex', 'fontsize', 16)
title('Residuals')
printFigure, print(fig, 'cylinder_Re100_res', '-dpdf', '-r0');
