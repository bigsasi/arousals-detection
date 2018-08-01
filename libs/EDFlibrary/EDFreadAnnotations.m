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
function [annotations, header] = EDFreadAnnotations(filename, doSkipTimeKeep)

annotations = [];

% Skips 'Time-keep' annotations
if isempty(doSkipTimeKeep)
    doSkipTimeKeep = 0;
end

general_header_size = 256; %bytes
one_signal_header_size = 256; %bytes

header = EDFreadHeader(filename);

if isempty(header)
    fprintf(1, 'Unable to read file header %s', filename);
    return;
end

fid = fopen(filename);

if (fid == -1)
    disp('Error opening the file');
    return;
end

bytes_full_data_record = 2 * sum([header.signals_info.num_samples_datarecord]);
header_length = general_header_size + header.num_signals * one_signal_header_size;

% Set the pointer to the annotations signal
signalIndex = header.num_signals;

if not(strcmpi(header.signals_info(signalIndex).label', 'EDF Annotations '))
    disp('Error: Not annotation signals was found in the file');
    return;
end

% Positioning on the begining of the selected file
fseek(fid, header_length + header.signals_info(signalIndex).signalOffset, 'bof');

record_size = header.signals_info(signalIndex).num_samples_datarecord;

% We will read all the annotations
numRecords2read = header.num_data_records;

% Data reading
data = fread(fid, record_size*numRecords2read, strcat(num2str(record_size),'*uint16=>uint16'), bytes_full_data_record - 2*record_size);
fclose(fid);

eventsall = typecast(data, 'uint8'); % From two's complement to 1-byte US-ASCII characters

tmp = find(eventsall == 20); % ASCII '20' separates each annotation after the time stamp (onset+duration)
% Detection of 'end of TAL': ASCII '20' followed by ASCII '0'. Each TAL may contain several annotations, and the last one finishes with this combination
% Note: Each TAL shares time_stamp and duration. Several TALS can be
%       included in the same block, but the first one contains an empty
%       annotation that indicates the time_stamp of the block itself.
numberTALends = tmp(eventsall(tmp+1) == 0); 

annotations.duration = [];
annotations.label = [];
annotations.offset = [];
for k = 1:length(numberTALends)
      
    if (k == 1)
        tal = eventsall(1:numberTALends(1));
    else
        tal = eventsall(numberTALends(k-1)+2:numberTALends(k));
    end
    if (ne(char(tal(1)),'+') && ne(char(tal(1)),'-'))
      % Remove previous characters til + or - (possibly by block transition, thus zeros at left are expected)
      strIdx = strfind(char(tal)', '+');
      if isempty(strIdx)
        strIdx = strfind(char(tal)', '-');
      end
      if isempty(strIdx)
        error('Incorrect start of TAL');
      end
      tal = tal(strIdx:end);
    end
   % Find offset and duration
   sep = find(tal == 20);
   indxDur = find(tal == 21);
   if isempty(indxDur)
       dur = 0;
       % TODO: the char() conversion is only consistent for the first 127
       % US-ASCII symbols (the ones mostly used), but starting 128 on, they
       % might be machine-dependent. For these the use of native2unicode
       % or equivalent function, specifying the UTF-8 encoding might be preferable
       % For the moment we leave it like this, as native2unicode is not
       % standarly supported by Octave
       % This note applies for the rest of similar conversions on this code
       %offset = str2double(native2unicode(tal(1:sep(1)-1))'); 
       offset = str2double(char(tal(1:sep(1)-1))');         
   else
       dur = str2double(char(tal(indxDur(1)+1:sep(1)-1))');
       offset = str2double(char(tal(1:indxDur(1)-1))');  
   end
   % Read annotations
   for k1 = 1:length(sep)-1
       tmp = char(tal(sep(k1)+1:sep(k1+1)-1));
 
       if isempty(tmp)
           annotations.label = [annotations.label, {'Time-keep'}];
       else
           %annotations.label = [annotations.label, textscan(tmp, '%s')];
           annotations.label = [annotations.label, {tmp'}];
       end
       annotations.offset = [annotations.offset, offset];
       annotations.duration = [annotations.duration, dur];
   end
end

if doSkipTimeKeep
    indxTimeKeep = strcmpi(annotations.label, 'Time-keep');
    annotations.offset(indxTimeKeep) = [];
    annotations.duration(indxTimeKeep) = [];
    annotations.label(indxTimeKeep) = [];
end