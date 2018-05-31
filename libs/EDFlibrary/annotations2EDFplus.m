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

function statusok = annotations2EDFplus(annotations, patId, recId, startdate, starttime, filename)

edf_block_size_limit = 61440; % In bytes

% Compute how many characters we will need to store
writenChars = 0;
num_necessary_blocks = 1; % At least one block is necessary

% Sort annotations before writing to file
[~,sortIndx] = sort([annotations.offset]);
annotations.offset = annotations.offset(sortIndx);
annotations.duration = annotations.duration(sortIndx);
annotations.label = annotations.label(sortIndx);

% Note: The time-keeping TAL for the first datablock it always assumes a zero onset, 
%       thus the corresponding TAL is: +0'20''20''0'
data = sprintf('%s%s', '+0', [20,20,0]);
writenChars = length(data);

% Each annotation will be included in a separate TAL
for k = 1:length(annotations.offset)
    
    dataForThisAnnotation = [];
    
    % Note: Each offset starts with '+' and finishes with unprintable ASCII '21'
    if (annotations.offset(k) >= 0)
      dataForThisAnnotation = [dataForThisAnnotation, sprintf('%s%s%s', '+', num2str(annotations.offset(k)), 21)];
    else
      dataForThisAnnotation = [dataForThisAnnotation, sprintf('%s%s', num2str(annotations.offset(k)), 21)];
    end
    % Note: Each duration finishes with unprintable ASCII '20'
    dataForThisAnnotation = [dataForThisAnnotation, sprintf('%s%s', num2str(annotations.duration(k)), 20)];
    % Note: Each annotation finishes with unprintable ASCII '20' and must not contain
    %       any '20' within. Also because we assume one TAL per annotation
    %       then the unprintable ASCII '0' is added to close the TAL
    dataForThisAnnotation = [dataForThisAnnotation, sprintf('%s%s', char(annotations.label(k)), [20,0])];
        
    % Notice each char "weights" 1-byte
    if (writenChars + length(dataForThisAnnotation) >= edf_block_size_limit)
        % This annotation needs to be set on a new block
        num_necessary_blocks = num_necessary_blocks + 1;
        
        % Current data is filled with zeros to the end of the datablock
        data = [data, sprintf('%s', zeros(1, edf_block_size_limit - writenChars))];
        
        % Reset the "writtenChars" counter
        writenChars = 0;
        
        % A a new time-keeping TAL is inserted before the annotation
        time_keep_tal = sprintf('%s%s%s', '+', num2str(annotations.offset(k)), [20,20,0]);
        dataForThisAnnotation = [time_keep_tal, dataForThisAnnotation]; 
    end
        
    data = [data, dataForThisAnnotation];
    writenChars = writenChars + length(dataForThisAnnotation);
        
    % TODO: An alternative might be to have instead "n" data pieces, thus data defined as
    % cell structure of "n" rows or similar, and then perform writing "per block"
    % and for each block just "write" the corresponding data{k}?
end

%Note: Unused bytes of the 'EDF Annotations' signal in the remainder of the data record are also filled with 0-bytes
data = [data, zeros(1, edf_block_size_limit - writenChars)];

fid = fopen(filename, 'w', 'ieee-le');

if (fid == -1)
    disp('Error creating output EDF+ file');
    return;
end

general_header_size = 256; %bytes
one_signal_header_size = 256; %bytes

% Write edf

% FIXED HEADER
header.version = 0;
header.local_patient_identification = patId;
header.local_recording_identification = recId;
header.startdate_recording = startdate;
header.starttime_recording = starttime;
header.num_signals = 1;
header.num_bytes_header = general_header_size + one_signal_header_size;
header.reserved = 'EDF+C';
header.duration_data_record = 0; % Note we assume an 'Annotations only' EDF+ file
% Each character "weights" 1-byte, thus the total number of necessary bytes
% equal the total number of characters to write
header.num_data_records = num_necessary_blocks;

fprintf(fid, trimAndFillWithBlanks(num2str(header.version), 8));   % version
fprintf(fid, '%-80s', header.local_patient_identification);
fprintf(fid, '%-80s', header.local_recording_identification);
fprintf(fid, '%-8s', header.startdate_recording);
fprintf(fid, '%-8s', header.starttime_recording);
actualHeaderBytes = general_header_size + one_signal_header_size*header.num_signals;
if (ne(actualHeaderBytes, header.num_bytes_header))
    disp('EDFwriteSignal: Warning, num_bytes_header does not match the actual number of header bytes. Fixed!');
end
fprintf(fid, trimAndFillWithBlanks(num2str(actualHeaderBytes), 8));
fprintf(fid, '%-44s', header.reserved);
fprintf(fid, trimAndFillWithBlanks(num2str(header.num_data_records), 8));
fprintf(fid, trimAndFillWithBlanks(num2str(header.duration_data_record), 8));
fprintf(fid, trimAndFillWithBlanks(num2str(header.num_signals), 4));

% SIGNAL DEPENDENT HEADER
header.signals_info(1).label = 'EDF Annotations';
header.signals_info(1).transducer_type = '';
header.signals_info(1).physical_dimension = '';
header.signals_info(1).physical_min = -32768;
header.signals_info(1).physical_max = 32767;
header.signals_info(1).digital_min = -32768;
header.signals_info(1).digital_max = 32767;
header.signals_info(1).prefiltering = '';
% One sample is a 2-byte (samples are encoded in EDF(+) using 16 bit pieces), and each character is using 1-byte, thus each character "weights" half a sample
header.signals_info(1).num_samples_datarecord = edf_block_size_limit/2;
header.signals_info(1).reserved = '';

fprintf(fid, '%-16s', header.signals_info(1).label);
fprintf(fid, '%-80s', header.signals_info(1).transducer_type);
fprintf(fid, '%-8s', header.signals_info(1).physical_dimension);
fprintf(fid, trimAndFillWithBlanks(num2str(header.signals_info(1).physical_min), 8));
fprintf(fid, trimAndFillWithBlanks(num2str(header.signals_info(1).physical_max), 8));
fprintf(fid, trimAndFillWithBlanks(num2str(header.signals_info(1).digital_min), 8));
fprintf(fid, trimAndFillWithBlanks(num2str(header.signals_info(1).digital_max), 8));
fprintf(fid, '%-80s', header.signals_info(1).prefiltering);
fprintf(fid, trimAndFillWithBlanks(num2str(header.signals_info(1).num_samples_datarecord), 8));
fprintf(fid, '%-32s', header.signals_info(1).reserved);

% DATA WRITING
header_length = general_header_size + header.num_signals * one_signal_header_size;
current_position = ftell(fid); % in bytes

if ne(header_length, current_position)
    disp('something wrong could be happening');
end

fwrite(fid, data, 'char');
statusok = (fclose(fid) == 0);
