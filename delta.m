function d = delta(r, dr)

    d = zeros(size(r));

    % Mask accounting for the discrete delta function in matrix form
    m3 = r > 1.5*dr;
    m1 = r <= 0.5*dr;    
    d1 = (1 + sqrt(1 - 3*(r.*m1./dr).^2))./(3*dr);
    m2 = ~(m1 | m3);
    d2 = (5 - 3*abs(r.*m2)./dr - sqrt(1 - 3*(1 - abs(r.*m2)./dr).^2))./(6*dr);
    
    d(m1) = d1(m1);
    d(m2) = d2(m2);
    
end

