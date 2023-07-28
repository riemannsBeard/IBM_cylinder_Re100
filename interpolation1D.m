function E_ = interpolation1D(xi, x, dx, normalized, nelem)

if nargin < 5
    nelem = 2;
end

[~, xi_x] = min(abs(xi - x), [], 2);

xi_dx = dx(xi_x);

if normalized == 1
    sf = xi_dx;
else
    sf = ones(size(xi_dx));
end

e_i = [];
e_j = [];
e_val = [];

for j = 1:length(xi)
    xj = xi_x(j);
    
    e_i = [e_i, j*ones(1, 2*nelem + 1)];
    e_j = [e_j, (xj - nelem):(xj + nelem)];
    deltaj = delta((xi(j) - x(xj - nelem:xj + nelem)), xi_dx(j));
    e_val = [e_val, sf(j)*deltaj];
end

e_i = e_i(:)';
e_j = e_j(:)';
e_val = e_val(:)';

E_ = sparse(e_i, e_j, e_val, length(xi), length(x));

end

