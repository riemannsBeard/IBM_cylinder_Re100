function d = delta(r, dr)

d = zeros(size(r));
absr = abs(r./dr);
m3 = absr > 1.5;
m1 = absr <= 0.5;

d1 = (1 + sqrt(1 - 3*m1.*absr.^2))./(3*dr);
m2 = ~(m1 | m3);
d2 = (5 - 3*absr.*m2 - sqrt(1 - 3*(1 - absr.*m2).^2))./(6*dr);

d = d1.*m1 + d2.*m2;

end

