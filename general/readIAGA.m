function data = readIAGA(filepath,varargin)
% readIAGA(filepath) loads the specified IAGA-2002 format file into a
% MATLAB structure
%
% See https://www.ngdc.noaa.gov/IAGA/vdat/IAGA2002/iaga2002format.html
% for an in-depth description of the IAGA 2002 data format
%
% data = readIAGA(filepath,'Format','XYZ')
%
%   Inputs:
%       -filepath: Path to file (string). Can be relative to current
%                  directory or a full path e.g.
%                  'data/cnb20220822qmin.min' or
%                  'C:/data/cnb20220822qmin.min'
%       -'Format': Optional argument. Type 'Format' followed by
%                  'XYZ', or 'HDZ'. Will force modify the IAGA-2002 file to
%                  the specified format. By default it will be in the
%                  format of the IAGA file.
%
%   Outputs:
%       -data: A MATLAB structure containing the data in neatly dot-indexed
%              fields e.g.
%         data = 
%            struct with fields:
%                         format: 'IAGA-2002'
%                 source_of_data: 'Indian Institute of Geomagnetism'
%                   station_name: 'Alibag'
%                      iaga_code: 'ABG'
%              geodetic_latitude: 18.6380
%             geodetic_longitude: 72.8720
%                      elevation: 7
%                       reported: 'XYZF'
%             sensor_orientation: 'HDZF'
%               digital_sampling: '1-second'
%             data_interval_type: '1-minute'
%                      data_type: 'variation'
%               publication_date: ''
%                       comments: '# This data file was created using...'
%                           data: [1439×7 table]
% 
% SEE ALSO: xyz2hdz, hdz2xyz
%
% Dr. Matthew Gard, 2022

% Default values and optional arguments
%---------------------------------------
% Default values
forceXYZ = 0;
forceHDZ = 0;
iagaMissing = 99999; % 99999 - no data, 88888 - not recorded
iagaNotRecorded = 88888;

% Parse optional arguments
for i = 1:length(varargin)
    % Check if format override given
    if strcmpi(varargin{i},'Format')
        if i < length(varargin)
            if strcmpi(varargin{i+1},'HDZ')
                forceHDZ = 1;
            elseif strcmpi(varargin{i+1},'XYZ')
                forceXYZ = 1;
            end
        end
    end
end


% Validation steps and initialising variables
%---------------------------------------------
%   filename must be a string
if ~isa(filepath, 'char') && ~isa(filepath, 'string')
    error(sprintf('File path must be a character array or string\nmagdata = iagaread(''path/to/file.format'')'));
end
%   filename must exist
if ~isfile(filepath);error('File ''%s'' not found.',filepath);end

% Initialise output structure
data = struct();
data.format = '';
data.source_of_data = '';
data.station_name = '';
data.iaga_code = '';
data.geodetic_latitude = NaN;
data.geodetic_longitude = NaN;
data.elevation = NaN;
data.reported = '';
data.sensor_orientation = '';
data.digital_sampling = '';
data.data_interval_type = '';
data.data_type = '';
data.publication_date = '';
data.comments = '';
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
[fID,message] = fopen(filepath,'rt');
if fID == -1
    error('Opening file failed: %s',message)
end

% Loop through file gathering header information and header line count
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
        if isempty(data.comments)
            data.comments = thisline;
        else
            data.comments = sprintf('%s\n%s',data.comments,thisline);
        end
    end

    % Loop through fields
    for i = 1:length(meta_data)
        if contains(thisline,meta_data{i})
            if meta_data_isnum(i)
                data.(meta_data_labels{i}) = str2double(...
                    strtrim(erase(thisline,meta_data{i})));
            else
                data.(meta_data_labels{i}) = strtrim(...
                    erase(thisline,meta_data{i}));
            end
        end
    end

    % Update header line count
    headerLines = headerLines + 1;
end

% Close file
fclose(fID);

% Grab the data below the headers
geomag_data = readtable(filepath,'FileType','text',...
    'HeaderLines',headerLines,'Format','%{uuuu-MM-dd}D%{HH:mm:ss.SSS}D%d%f%f%f%f');

% Grab column names from reported field e.g. XYZF vs HDZF
if ~isempty(data.reported)
    columndatalabels = repmat({''}, 1, length(data.reported));
    for i = 1:length(data.reported)
        columndatalabels{i} = upper(data.reported(i));
    end
    columndatalabels = [{'date'},{'time'},{'DOY'},columndatalabels];
    geomag_data.Properties.VariableNames = columndatalabels;
    geomag_data.time = timeofday(geomag_data.time);
else
    % Try and pull from column labels if reported format is missing
    % But it shouldnt be
    % thisline is the last line before data, should still be in memory
    if contains(regexprep(thisline,'\s+',' '),'DATE TIME DOY')
        thisline = strsplit(thisline);
        % Remove the DATE, TIME and DOY
        thisline = thisline(4:end);
    end
    columndatalabels = repmat({''}, 1, length(thisline));
    for i = 1:length(thisline)
        % Format of colum names should be SSSX, SSSY, SSSZ
        columndatalabels{i} = upper(thisline{i}(4));
    end
    columndatalabels = [{'date'},{'time'},{'DOY'},columndatalabels];
    geomag_data.Properties.VariableNames = columndatalabels;
    geomag_data.time = timeofday(geomag_data.time);
end

% If format wanted is HDZ/XYZ
if forceXYZ
    % Requires H and D
    if any(ismember(geomag_data.Properties.VariableNames,'H')) && any(ismember(geomag_data.Properties.VariableNames,'D'))
        % Find the indices where the data is missing and preserve this
        indMissing = geomag_data.H==iagaMissing | geomag_data.D==iagaMissing;
        indNotRecoded = geomag_data.H==iagaNotRecorded | geomag_data.D==iagaNotRecorded;

        [geomag_data.H,geomag_data.D] = hdz2xyz(geomag_data.H,geomag_data.D);
        idx1 = find(ismember(geomag_data.Properties.VariableNames,'H'));
        idx2 = find(ismember(geomag_data.Properties.VariableNames,'D'));
        geomag_data.Properties.VariableNames([idx1,idx2]) = {'X','Y'};

        % Restore missing data values
        geomag_data{indMissing,{'X','Y'}} = iagaMissing;
        geomag_data{indNotRecoded,{'X','Y'}} = iagaNotRecorded;

    elseif any(ismember(geomag_data.Properties.VariableNames,'X')) || any(ismember(geomag_data.Properties.VariableNames,'Y'))
        % Do nothing and exit
    else
        fprintf('Unable to force XYZ format for file %s (Missing H or D data, or format incorrect)\n',filepath)
    end
end
if forceHDZ
    % Requires X and Y
    if any(ismember(geomag_data.Properties.VariableNames,'X')) && any(ismember(geomag_data.Properties.VariableNames,'Y'))
        % Find the indices where the data is missing and preserve this
        indMissing = geomag_data.X==iagaMissing | geomag_data.Y==iagaMissing;
        indNotRecoded = geomag_data.X==iagaNotRecorded | geomag_data.Y==iagaNotRecorded;

        [geomag_data.X,geomag_data.Y] = xyz2hdz(geomag_data.X,geomag_data.Y);
        idx1 = find(ismember(geomag_data.Properties.VariableNames,'X'));
        idx2 = find(ismember(geomag_data.Properties.VariableNames,'Y'));
        geomag_data.Properties.VariableNames([idx1,idx2]) = {'H','D'};

        % Restore missing data values
        geomag_data{indMissing,{'H','D'}} = iagaMissing;
        geomag_data{indNotRecoded,{'H','D'}} = iagaNotRecorded;
    elseif any(ismember(geomag_data.Properties.VariableNames,'H')) || any(ismember(geomag_data.Properties.VariableNames,'D'))
        % Do nothing and exit
    else
        fprintf('Unable to force HDZ format for file %s (Missing X or Y data, or format incorrect)\n',filepath)
    end
end

% Save the geomag data to .data field
data.data = geomag_data;

return