function I = inclination(H,Z)
% inclination(H,Z) calculates magnetic inclination given H and Z magnetic
% field values according to the function: Inclination I = tan-1(Z/H)
% 
%   I = inclination(H,Z)
%
%   Inputs:
%       -H: Horizontal intensity (H) of the magnetic field
%       -Z: Vertical intensity (Z) of the magnetic field
%   H and Z should be arrays of numeric data of equal size
% 
%   Outputs:
%       -I: Inclination (I) of the magnetic field (in degrees)
%   I will be the same size as the inputs
%
% SEE ALSO: xyz2hdz
% 
% Dr. Matthew Gard, 2022

% Validate inputs
%-----------------
if (~isnumeric(H))
    error('Input variable H is not numeric')
end
if (~isnumeric(Z))
    error('Input variable Z is not numeric')
end

if ~isequal(size(H),size(Z))
    error('Input variables H and Z are not the same size')
end

% Processing
%-----------------
I = atand(Z./H);

return