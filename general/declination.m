function D = declination(X,Y)
% declination(X,Y) calculates magnetic declination given X and Y magnetic
% field values according to the function: Declination D = tan-1(Y/X)
% 
%   D = declination(X,Y)
%
%   Inputs:
%       -X: North component (X) of the magnetic field
%       -Y: East component (Y) of the magnetic field
%   X and Y should be arrays of numeric data of equal size
% 
%   Outputs:
%       -D: Declination (D) of the magnetic field (in degrees)
%   D will be the same size as the inputs
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
D = atand(Y./X);

return