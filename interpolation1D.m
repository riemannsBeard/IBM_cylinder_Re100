function E_ = interpolation1D(xi, x, dx, normalized)

    [~, xi_x] = min(abs(xi - x)); 
    xi_dx = dx(xi_x);

    if normalized == 1
        sf = xi_dx;
    else
        sf = ones(size(xi_dx));
    end
        
    
end

