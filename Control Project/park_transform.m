function [vd, vq] = park_transform(va, vb, theta)
% Transformada de Park - va/vb para vd/vq
vd = va * cos(theta) + vb * sin(theta);
vq = -va * sin(theta) + vb * cos(theta);
end