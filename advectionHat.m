function [ Nhat ] = advectionHat(grid, u, v, Nx, Ny, bc)
   
    u = reshape(u, Ny, Nx-1);
    v = reshape(v, Ny-1, Nx);
    
    %% Condiciones de contorno (en los nodos de la malla externa)
    
%     % u
%     bc.uN = grid.X*0;           
%     bc.uS = grid.X*0;               
%     bc.uW = 0.5*(grid.Y(2:end) + grid.Y(1:end-1))'*0 + 1;
%     bc.uE = 0.5*(grid.Y(2:end) + grid.Y(1:end-1))'*0 + 1;
%     
%     % v
%     bc.vN = 0.5*(grid.X(2:end) + grid.X(1:end-1))*0;
%     bc.vS = 0.5*(grid.X(2:end) + grid.X(1:end-1))*0;
%     bc.vW = grid.Y'*0;
%     bc.vE = grid.Y'*0;
    
    %% Tratamiento de los términos convectivos
    
    % Añado contornos laterales y filas fantasma a la malla interna de las
    % velocidades
%     ue = [bc.uW u bc.uE];
%     ue = [2*bc.uS - ue(1,:); ue; 2*bc.uN - ue(end,:)];
%     
%     ve = [bc.vS; v; bc.vN];
%     ve = [2*bc.vW - ve(:,1) ve 2*bc.vE - ve(:,end)];
% 
%     % Interpolacion lineal por filas para que ambas velocidades esten en
%     % los nodos de la malla externa
%     ua = 0.5*(ue(2:end,:) + ue(1:end-1,:));
%     va = 0.5*(ve(:,2:end) + ve(:,1:end-1));
    
%     ud = diff(ue)/2;
%     vd = diff(ve')'/2;
    
%     % Calculo las derivadas que involucran al producto de ambas velocidades
%     Nhat.vux = diff(ua.*va, 1, 2)./grid.dX';
%     Nhat.uvy = diff(ua.*va, 1, 1)./grid.dY;
% 
%     % Interpolacion lineal por columnas para que, al derivar, ambas
%     % velocidades queden en la posicion de la presion
%     ua = 0.5*(ue(2:end-1,1:end-1) + ue(2:end-1, 2:end));
%     va = 0.5*(ve(1:end-1,2:end-1) + ve(2:end, 2:end-1));

    Nhat.u = zeros(size(u));
    Nhat.v = zeros(size(v));
    
    % Calculo de las derivadas que involucran al termino de la velocidad al
    % cuadrado
    u2 = u.*u;
    v2 = v.*v;
    
    u2a = zeros(Ny, Nx);
    v2a = zeros(Ny, Nx);
    
    u2a(:,1) = 0.5*(bc.uW.^2 + u2(:,1));
    u2a(:,2:end-1) = 0.5*(u2(:,2:end) + u2(:,1:end-1));
    u2a(:,end) = 0.5*(bc.uE.^2 + u2(:,end));
    
    v2a(1,:) = 0.5*(bc.vS.^2 + v2(1,:));
    v2a(2:end-1,:) = 0.5*(v2(2:end,:) + v2(1:end-1,:));
    v2a(end,:) = 0.5*(bc.vN.^2 + v2(end,:));
    
    Nhat.u2x = diff(u2a, 1, 2)./grid.dXp';
    Nhat.v2y = diff(v2a, 1, 1)./grid.dYp;

    Nhat.uv = 0.25*(u(2:end,:) + u(end,:)).*(v(:,2:end) + v(:,end-1));
    
    Nhat.uvS = 0.5*bc.uS.*(bc.vS(:,2:end) + bc.vS(:,1:end-1));
    Nhat.uvN = 0.5*bc.uN.*(bc.vN(:,2:end,:) + bc.vN(:,1:end-1));
    
    Nhat.uvW = 0.5*bc.vW.*(bc.uW(2:end,:) + bc.uW(1:end-1,:));
    Nhat.uvE = 0.5*bc.vE.*(bc.uE(2:end,:) + bc.uE(1:end-1,:));    
    
    %% Vectorizacion y normalizacion
    
    Nhat.u = Nhat.u + diff([Nhat.uvS; Nhat.uv; Nhat.uvN], 1, 1)./grid.dY;
    Nhat.v = Nhat.v + diff([Nhat.uvW Nhat.uv Nhat.uvE], 1, 2)./grid.dX';

    Nhat.u = reshape(Nhat.u, [], 1);
    Nhat.v = reshape(Nhat.v, [], 1);
      
    %% Velocidades en la malla externa
    
%     ua = (ue(2:end,:) + ue(1:end-1,:))/2;
%     va = (ve(:,2:end) + ve(:,1:end-1))/2;
    
end

