function bwr = bluewhitered(m)
% bluewhitered
%   bluewhitered(M) returns an M-by-3 matrix containing the red blue colormap
%   The colors begin with blue, transitioning to white, and then to red.
%   bluewhitered, by itself, is the same length as the current figure's colormap.
%   If no figure exists, MATLAB uses the length of the default colormap.
%
% Dr. Matthew Gard, 2022

if nargin < 1
   f = get(groot,'CurrentFigure');
   if isempty(f)
      m = size(get(groot,'DefaultFigureColormap'),1);
   else
      m = size(f.Colormap,1);
   end
end

% Force m positive integer
m = ceil(abs(m));

% Solid blue:    0 0 1
% Solid white:   1 1 1
% Solid red:     1 0 0

% Must be greater than or equal to 3, otherwise just set it to rwb
% If you're using rwb you want at least 3 colors!
if m > 3
    % If m is even, double the size of the white portion i.e. end bw at
    % [1,1,1] and also start wr at [1,1,1]
    if mod(m,2)==0 % iseven
        bw = [linspace(0,1,m/2)', linspace(0,1,m/2)',...
            ones(m/2,1)];
        wr = [ones(m/2,1), linspace(1,0,m/2)',...
            linspace(1,0,m/2)'];
        bwr = [bw;wr];
    else
        % If m is not even, do each division piece wise without the
        % duplication of the central white space
        lowerm = floor(m/2);
        upperm = ceil(m/2);
        bw = [linspace(0,1-(1/lowerm),lowerm)', linspace(0,1-(1/lowerm),lowerm)',...
            ones(lowerm,1)];
        wr = [ones(upperm,1), linspace(1,0,upperm)',...
            linspace(1,0,upperm)'];
        bwr = [bw;wr];
    end
else
    bwr = [0 0 1;1 1 1;1 0 0];
end

return