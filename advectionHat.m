function [ Nhat, ua, va ] = advectionHat(grid, u, v, Nx, Ny, bc)
   
    u = reshape(u, Ny, Nx-1);
    v = reshape(v, Ny-1, Nx);
          

    %% Tratamiento de los términos convectivos
    
    Nhat.u = zeros(size(u));
    Nhat.v = zeros(size(v));
    
    % Calculo de las derivadas que involucran al termino de la velocidad al
    % cuadrado
    u2 = u.*u;
    v2 = v.*v;
    
    u2c = zeros(Ny, Nx);
    v2c = zeros(Ny, Nx);
    
    u2c(:,1) = 0.5*(bc.uW.^2 + u2(:,1));
    u2c(:,2:end-1) = 0.5*(u2(:,2:end) + u2(:,1:end-1));
    u2c(:,end) = 0.5*(bc.uE.^2 + u2(:,end));
        
    Nhat.u = diff(u2c')'./grid.dXp';
    
    v2c(1,:) = 0.5*(bc.vS.^2 + v2(1,:));
    v2c(2:end-1,:) = 0.5*(v2(2:end,:) + v2(1:end-1,:));
    v2c(end,:) = 0.5*(bc.vN.^2 + v2(end,:));
    
    Nhat.v = diff(v2c)./grid.dYp;

    % uv
    uv = 0.25*(u(2:end,:) + u(1:end-1,:)).*(v(:,2:end) + v(:,1:end-1));
    
       
    uvS = 0.5*bc.uS.*(bc.vS(2:end) + bc.vS(1:end-1));
    uvN = 0.5*bc.uN.*(bc.vN(2:end) + bc.vN(1:end-1));
    
    uvW = 0.5*bc.vW.*(bc.uW(2:end) + bc.uW(1:end-1));
    uvE = 0.5*bc.vE.*(bc.uE(2:end) + bc.uE(1:end-1));   
    
    %% Vectorizacion y normalizacion
    
    Nhat.u = Nhat.u + diff([uvS; uv; uvN], 1, 1)./grid.dY;
    Nhat.v = Nhat.v + diff([uvW, uv, uvE], 1, 2)./grid.dX';
    
    Nhat.u = reshape(Nhat.u, [], 1);
    Nhat.v = reshape(Nhat.v, [], 1);

%%    
      
    % u
    bc.uN = grid.X*0;           
    bc.uS = grid.X*0;               
    bc.uW = 0.5*(grid.Y(2:end) + grid.Y(1:end-1))'*0 + 1;
    bc.uE = 0.5*(grid.Y(2:end) + grid.Y(1:end-1))'*0 + 1;
    
    % v
    bc.vN = 0.5*(grid.X(2:end) + grid.X(1:end-1))*0;
    bc.vS = 0.5*(grid.X(2:end) + grid.X(1:end-1))*0;
    bc.vW = grid.Y'*0;
    bc.vE = grid.Y'*0;
    
    %% Tratamiento de los términos convectivos
    
    % Añado contornos laterales y filas fantasma a la malla interna de las
    % velocidades
    ue = [bc.uW u bc.uE];
    ue = [2*bc.uS - ue(1,:); ue; 2*bc.uN - ue(end,:)];
    
    ve = [bc.vS; v; bc.vN];
    ve = [2*bc.vW - ve(:,1) ve 2*bc.vE - ve(:,end)];

    % Interpolacion lineal por filas para que ambas velocidades esten en
    % los nodos de la malla externa
    ua = 0.5*(ue(2:end,:) + ue(1:end-1,:));
    va = 0.5*(ve(:,2:end) + ve(:,1:end-1));
        
end

