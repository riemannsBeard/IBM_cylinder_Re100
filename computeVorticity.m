function vorticity = computeVorticity(ua, va, grid)

vorticity = gradient(va)./gradient(grid.x) - gradient(ua')'./gradient(grid.y')';

end

