function E_ = interpolation1D(xi, x, dx, normalized, nelem)

    [~, xi_x] = min(abs(xi - x)); 
    xi_dx = dx(xi_x);
    
    e_i = [];
    e_j = [];
    e_val = [];

    if normalized == 1
        sf = xi_dx;
    else
        sf = ones(size(xi_dx));
    end
    
    for j = 1:length(xi)
        e_i = [e_i; j*ones(2*nelem + 1)];
        e_j = [e_j; xj-nelem:xj+nelem+1];
        deltaj = delta((xi(j) - x(xj-nelem:xj+nelem+1)), xi_dx(j));
        e_val = [e_val; sf(j)*deltaj];
    end
    
    e_i = e_i(:);
    e_j = e_j(:);
    e_val = e_val(:);
    
    E_ = sparse(e_val, [e_i, e_j], length(xi), length(x));
        
end

