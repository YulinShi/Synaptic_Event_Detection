function somanew = somaPositionTransformer_Peak_Detection(somaX, somaY, spatialRotation, xPatternOffset, yPatternOffset)
% somaPositionTransformer

% editing:
% gs may 2006
% --------------------------------------------

if nargin == 0
    somaX = -9;
    somaY = 223;
    spatialRotation = -4;
    xPatternOffset = 0;
    yPatternOffset = -100;
end
somaXoffset = somaX - xPatternOffset;
somaYoffset = somaY - yPatternOffset;
rotationAngleRadians = (-1) * spatialRotation * (pi/180);
[theta, rho] = cart2pol(somaXoffset, somaYoffset);
[somaXnew, somaYnew] = pol2cart(theta + rotationAngleRadians, rho);
somanew=[somaXnew, somaYnew];
% disp([num2str(round(somaXnew)) '    ' num2str(round(somaYnew))]);