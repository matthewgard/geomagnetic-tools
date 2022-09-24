function [H,D] = xyz2hdz(X,Y)
% xyz2hdz(X,Y) computes the Horizontal Intensity (H) and Declination (D)
% components of the magnetic field given xy data. Z is not necessary for
% the calculation.
%   Horizontal  H = sqrt(x2+y2)
%   Declination D = tan-1(Y/X)
% 
%   [H,D] = xyz2hdz(X,Y)
%
%   Inputs:
%       -X: North component (X) of the magnetic field
%       -Y: East component (Y) of the magnetic field
%   X and Y should be arrays of numeric data of equal size
% 
%   Outputs:
%       -H: Horizontal intensity (H) of the magnetic field
%       -D: Declination (D) of the magnetic field (in degrees)
%   H and D will be the same size as the inputs X and Y
%
% SEE ALSO: hdz2xyz
%
% Dr. Matthew Gard, 2022

% Validate inputs
%-----------------
if (~isnumeric(X))
    error('Input variable X is not numeric')
end
if (~isnumeric(Y))
    error('Input variable Y is not numeric')
end

if ~isequal(size(X),size(Y))
    error('Input variables X and Y are not the same size')
end

% Processing
%-----------------
H = sqrt(X.^2+Y.^2);
D = atand(Y./X);

return