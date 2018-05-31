%% The code below is based on the methods described in the following reference(s):
% 
% [1] - I. Fernández-Varela, D. Alvarez-Estevez, E. Hernández-Pereira, V. Moret-Bonillo, 
% "A simple and robust method for the automatic scoring of EEG arousals in
% polysomnographic recordings", Computers in Biology and Medicine, vol. 87, pp. 77-86, 2017 
%
% Copyright (C) 2017 Isaac Fernández-Varela
% Copyright (C) 2017 Diego Alvarez-Estevez

%% This program is free software: you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation, either version 3 of the License, or
%% (at your option) any later version.

%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.

%% You should have received a copy of the GNU General Public License
%% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [edfSignal, header] = EDFreadSignal(filename, signalIndex, startRec, blocks2read)

general_header_size = 256; %bytes
one_signal_header_size = 256; %bytes

header = EDFreadHeader(filename);

if isempty(startRec)
    startRec = 0;
end

fid = fopen(filename);

if (fid == -1)
    disp('Error opening the file');
    return;
end

bytes_full_data_record = 2 * sum([header.signals_info.num_samples_datarecord]);
header_length = general_header_size + header.num_signals * one_signal_header_size;

% Positioning on the begining of the selected file
fseek(fid, header_length + (bytes_full_data_record * startRec) + header.signals_info(signalIndex).signalOffset, 'bof');

record_size = header.signals_info(signalIndex).num_samples_datarecord;

if isempty(blocks2read)
    % By default we read the entire EDF
    numRecords2read = header.num_data_records;
else
    numRecords2read = blocks2read;
end

% Data reading
data = fread(fid, record_size*numRecords2read, strcat(num2str(record_size),'*uint16=>uint16'), bytes_full_data_record - 2*record_size);

fclose(fid);

data = typecast(data, 'int16'); % From two's complement to signed integers

Pmin = header.signals_info(signalIndex).physical_min;
Pmax = header.signals_info(signalIndex).physical_max;
Dmin = header.signals_info(signalIndex).digital_min;
Dmax = header.signals_info(signalIndex).digital_max;

logConversion = strcmp(header.signals_info(signalIndex).physical_dimension', 'Filtered');
if logConversion
    % Parse prefiltering to set this values
    tokens = textscan(header.signals_info(signalIndex).prefiltering, 'sign*LN[sign*(at %.1fHz)/(%.5f)]/(%.5f)(Kemp:J Sleep Res 1998-supp2:132)');
    if isempty(tokens{2})
        disp('Warning: assigned default values to logConversion');
        LogFloatY0 = 0.0001; % default value
    else
        LogFloatY0 = tokens{2};
    end
    if isempty(tokens{3})
        disp('Warning: assigned default values to logConversion');
        LogFloatA = 0.001; % default value
    else
        LogFloatA = tokens{3};
    end
    edfSignal = Pmin + (Pmax - Pmin) * (ExpIntegerVector(double(data), LogFloatY0, LogFloatA) - Dmin)/(Dmax - Dmin);
else
    edfSignal = Pmin + (Pmax - Pmin) * (double(data) - Dmin)/(Dmax - Dmin);
end

