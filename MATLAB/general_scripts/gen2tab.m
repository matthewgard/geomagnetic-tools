function output_table = gen2tab(gen_file)
% gen2tab(gen_file) loads a DGGRID into a MATLAB table
%
% DGGRID is a command-line application for generating and manipulating
% icosahedral discrete global grids (DGGs).
% https://github.com/sahrk/DGGRID
%
% output_table = gen2tab(gen_file)
%
%   Inputs:
%       -gen_file: DGGRID generated .gen file path
%
%   Outputs:
%       -output_table: Loaded DGGRID file as a MATLAB table
%
% output_table table format: ind, coordinates, vertices
%  ind (integer):               Index
%  coordinates (1x2 vector):    Latitude and Longitude
%  vertices (Xx2 cell array):   Vertices surrounding the central coordinate
% 
% SEE ALSO: 
%
% Dr. Matthew Gard, 2022

% Validation steps and initialising variables
%---------------------------------------------
%   filename must be a string
if ~isa(gen_file, 'char') && ~isa(gen_file, 'string')
    error(sprintf('File path must be a character array or string\nmagdata = iagaread(''path/to/file.format'')'));
end
%   filename must exist
if ~isfile(gen_file);error('File ''%s'' not found.',gen_file);end

% Initialise output table
output_table.ind = [];
output_table.coordinates = [];
output_table.vertices = [];

% Open the file
fprintf('Loading DGGRID file %s (This may take some time)\n',gen_file)
fid = fopen(gen_file,'rt');
i = 1;
% Loop through lines
tline = fgetl(fid);
while ischar(tline)
    line_var = str2double(strsplit(tline));
    if length(line_var) > 2
        i = line_var(1);
        output_table(i).coordinates = [line_var(3) line_var(2)];
        output_table(i).ind = i;
        tline = fgetl(fid);
    elseif length(line_var) == 2
        output_table(i).vertices(end+1,:) = [line_var(2) line_var(1)];
        output_table(i).ind = i;
        tline = fgetl(fid);
    elseif isnan(line_var)
        tline = fgetl(fid);
        continue
    end
end

fprintf('Number of positions loaded: %i\n',length(output_table))
output_table = struct2table(output_table);

return