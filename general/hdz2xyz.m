function [X,Y] = hdz2xyz(H,D)
% hdz2xyz(H,D) computes the North (X) and East (Y) components of the 
% magnetic field given xy data. Z is not necessary for the calculation.
%   North X = Hcos(D)
%   East Y = Hsin(D)
% 
%   [X,Y] = hdz2xyz(H,D)
%
%   Inputs:
%       -H: Horizontal intensity (H) of the magnetic field
%       -D: Declination (D) of the magnetic field (in degrees)
%   H and D should be arrays of numeric data of equal size
% 
%   Outputs:
%       -X: North component (X) of the magnetic field
%       -Y: East component (Y) of the magnetic field
%   X and Y will be the same size as the inputs H and D
%
% SEE ALSO: xyz2hdz
%
% Dr. Matthew Gard, 2022

% Validate inputs
%-----------------
if (~isnumeric(H))
    error('Input variable H is not numeric')
end
if (~isnumeric(D))
    error('Input variable D is not numeric')
end

if ~isequal(size(H),size(D))
    error('Input variables H and D are not the same size')
end

% Processing
%-----------------
X = H.*cosd(D);
Y = H.*sind(D);

return