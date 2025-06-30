function cor = cen2cor(rcen,ccen,half_det)
% Convert from center to corner reference
rcor = half_det - rcen;
ccor = half_det + ccen;
cor = [rcor ccor];
end

