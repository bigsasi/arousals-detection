%% The code below is based on the methods described in the following reference(s):
% 
% [1] - I. Fernández-Varela, D. Alvarez-Estevez, E. Hernández-Pereira, V. Moret-Bonillo, 
% "A simple and robust method for the automatic scoring of EEG arousals in
% polysomnographic recordings", Computers in Biology and Medicine, vol. 87, pp. 77-86, 2017
%
% [2] - D. Alvarez-Estevez, I. Fernández-Varela, "Large-scale validation of an automatic EEG arousal detection
% algorithm using different heterogeneous databases", Sleep Medicine, vol. 57, pp. 6-14, 2019 
%
% Copyright (C) 2017-2019 Isaac Fernández-Varela
% Copyright (C) 2017-2019 Diego Alvarez-Estevez

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

function header = EDFreadHeader(filename)

fid = fopen(filename);

if (fid == -1)
    header = [];
    disp('Error opening the file');
    return;
end

general_header_size = 256; %bytes
one_signal_header_size = 256; %bytes

header.version = str2double(fread(fid, 8, 'uchar=>char')');
header.local_patient_identification = fread(fid, 80, 'uchar=>char');
header.local_recording_identification = fread(fid, 80, 'uchar=>char');
header.startdate_recording = fread(fid, 8, 'uchar=>char');
header.starttime_recording = fread(fid, 8, 'uchar=>char');
header.num_bytes_header = str2double(fread(fid, 8, 'uchar=>char')');
header.reserved = fread(fid, 44, 'uchar=>char');
header.num_data_records = str2double(fread(fid, 8, 'uchar=>char')');
header.duration_data_record = str2double(fread(fid, 8, 'uchar=>char')'); %in seconds
header.num_signals = str2double(fread(fid, 4, 'uchar=>char')');

ns = header.num_signals;

for k = 1:ns
    header.signals_info(k).label = fread(fid, 16, 'uchar=>char');
end
for k = 1:ns
    header.signals_info(k).transducer_type = fread(fid, 80, 'uchar=>char');
end
for k = 1:ns
    header.signals_info(k).physical_dimension = fread(fid, 8, 'uchar=>char');
end
for k = 1:ns
    header.signals_info(k).physical_min = str2double(fread(fid, 8, 'uchar=>char')');
end
for k = 1:ns
    header.signals_info(k).physical_max = str2double(fread(fid, 8, 'uchar=>char')');
end
for k = 1:ns
    header.signals_info(k).digital_min = str2double(fread(fid, 8, 'uchar=>char')');
end
for k = 1:ns
    header.signals_info(k).digital_max = str2double(fread(fid, 8, 'uchar=>char')');
end
for k = 1:ns
    header.signals_info(k).prefiltering = fread(fid, 80, 'uchar=>char');
end
signalOffsets = zeros(1, ns); % In bytes
for k = 1:ns
    header.signals_info(k).num_samples_datarecord = str2double(fread(fid, 8, 'uchar=>char')');
    % NOTE: The two following are not specific EDF header fields, but are practical for EDF handling
    if (header.duration_data_record > 0)
        header.signals_info(k).sample_rate = header.signals_info(k).num_samples_datarecord / header.duration_data_record;
    else
        header.signals_info(k).sample_rate = [];
    end
    if (k > 1)
        signalOffsets(k) = signalOffsets(k - 1) + 2 * header.signals_info(k - 1).num_samples_datarecord;
    end
    header.signals_info(k).signalOffset = signalOffsets(k);
end
for k = 1:ns
    header.signals_info(k).reserved = fread(fid, 32, 'uchar=>char');
end

header_length = general_header_size + ns * one_signal_header_size;
current_position = ftell(fid); % in bytes

if ne(header_length, current_position)
    disp('something wrong could be happening');
end

fclose(fid);