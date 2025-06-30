function cen = cor2cen(rcor,ccor,half_det)
% Convert from corner to center reference
rcen = half_det - rcor;
ccen = ccor - half_det;
cen = [rcen ccen];
end

