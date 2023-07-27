function d = delta(r, dr)

    d = zeros(size(r));
    absr = abs(r./dr);

    % Mask accounting for the discrete delta function in matrix form
    m3 = absr > 1.5;
    m1 = absr <= 0.5;
    m2 = ~(m1 || m3);
    
    d(m1) = (1 + sqrt(1 - 3*absr(m1).^2))./(3*dr);
    d(m2) = (5 - 3*absr(m2) - sqrt(1 - 3*(1 - absr(m2)).^2))./(6*dr);
    

end

