function magdata = readIAGA(filename,varargin)
% readIAGA loads an IAGA 2002 format file into a struct data type
% 
%   magdata = iagaread('path/to/file.format')
%
%   See https://www.ngdc.noaa.gov/IAGA/vdat/IAGA2002/iaga2002format.html
%   for an in-depth description of the IAGA 2002 data format
%
%   Inputs:
%       -filename:
% 
%   Outputs:
%       -magdata:
%
% Dr. Matthew Gard, 2022




% for i = 1:length(varargin)
%     if strcmpi(varargin{i},'XYZ')
%         forceXYZ = 1;
%         forceHDZ = 0;
%     elseif strcmpi(varargin{i},'HDZ')
%         forceHDZ = 1;
%         forceXYZ = 1;
%     end
% end



% Validation steps and initialising variables
%---------------------------------------------
%   filename must be a string
if ~isa(filename, 'char') && ~isa(filename, 'string')
    error(sprintf('File path must be a character array or string\nmagdata = iagaread(''path/to/file.format'')'));
end
%   filename must exist
if ~isfile(filename);error('File ''%s'' not found.',filename);end

% Initialise output structure
magdata = struct();
magdata.format = '';
magdata.source_of_data = '';
magdata.station_name = '';
magdata.iaga_code = '';
magdata.geodetic_latitude = NaN;
magdata.geodetic_longitude = NaN;
magdata.elevation = NaN;
magdata.reported = '';
magdata.sensor_orientation = '';
magdata.digital_sampling = '';
magdata.data_interval_type = '';
magdata.data_type = '';
magdata.publication_date = '';
magdata.comments = '';


meta_data = [{'Format'},{'Source of Data'},{'Station Name'},...
        {'IAGA CODE'},{'Geodetic Latitude'},{'Geodetic Longitude'},...
        {'Elevation'},{'Reported'},{'Sensor Orientation'},...
        {'Digital Sampling'},{'Data Interval Type'},{'Data Type'},...
        {'Publication Date'}];
meta_data_isnum = [0,0,0,...
    0,1,1,...
    1,0,0,...
    0,0,0,...
    0];
meta_data_labels = [{'format'},{'source_of_data'},{'station_name'},...
        {'iaga_code'},{'geodetic_latitude'},{'geodetic_longitude'},...
        {'elevation'},{'reported'},{'sensor_orientation'},...
        {'digital_sampling'},{'data_interval_type'},{'data_type'},...
        {'publication_date'}];

% Loading data file
%-------------------
% Open file
[fID,message] = fopen(filename,'rt');
if fID == -1
    error('Opening file failed: %s',message)
end

% Loop through file gathering header information
headerLines = 0;
while true
    % Get the current line
    thisline = fgetl(fID);
    if ~ischar(thisline)
        % End of file, exit
        break
    end

    % Remove vertical bar at end if exists and whitespace
    thisline(strfind(thisline,'|')) = '';
    thisline = strtrim(thisline);


    % Check if empty line
    if isempty(thisline)
        headerLines = headerLines + 1;
        continue
    end

    % Check if data fields starting
    if contains(regexprep(thisline,'\s+',' '),'DATE TIME DOY')
        headerLines = headerLines + 1;
        break
    end

    % Check if line is comment, if yes, append/add to comments
    if strcmp(thisline(1),'#')
        if isempty(magdata.comments)
            magdata.comments = thisline;
        else
            magdata.comments = sprintf('%s\n%s',magdata.comments,thisline);
        end
    end

    % Loop through fields
    for i = 1:length(meta_data)
        if contains(thisline,meta_data{i})
            if meta_data_isnum(i)
                magdata.(meta_data_labels{i}) = str2double(...
                    strtrim(erase(thisline,meta_data{i})));
            else
                magdata.(meta_data_labels{i}) = strtrim(...
                    erase(thisline,meta_data{i}));
            end
        end
    end

    % Update header line count
    headerLines = headerLines + 1;
end

% More efficiently grab the data below the headers
% Get column names
geomag_data = readtable(filename,'FileType','text',...
    'HeaderLines',headerLines+1,'Format','%{uuuu-MM-dd}D%{HH:mm:ss.SSS}D%d%f%f%f%f');

if ~isempty(magdata.reported)
    columndatalabels = repmat({''}, 1, length(magdata.reported));
    for i = 1:length(magdata.reported)
        columndatalabels{i} = upper(magdata.reported(i));
    end
    columndatalabels = [{'date'},{'time'},{'DOY'},columndatalabels];
    geomag_data.Properties.VariableNames = columndatalabels;
    geomag_data.time = timeofday(geomag_data.time);
    magdata.data = geomag_data;
    
else
    % Try and pull from column labels?
    error('Error for now - column labels missing i.e. reported')
end



% Close file
fclose(fID);

return