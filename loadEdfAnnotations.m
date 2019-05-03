%LOADANNOTATIONS load annotations from a edf file. 
%   filename    complete file name (with path) 
%   annotations struct with a field for each event type. Field names are the 
%               ones used in the xml files, without white spaces.

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

function [annotations] = loadEdfAnnotations(file, fileAnnotations)
    logger = logging.Logger.getLogger('LoadAnnotations');
    
    header = EDFreadHeader(file);
    [edfAnnotations, headerAnnotations] = EDFreadAnnotations(fileAnnotations, 0);

    edfStart = stringToSeconds(header.starttime_recording'); 
    annotationsStart = stringToSeconds(headerAnnotations.starttime_recording');

    deltaTime = edfStart - annotationsStart;
    
    if deltaTime ~= 0
        logger.info('Annotations displaced: %d seconds', deltaTime);
    end
    
%     annotations = [];
    annotations.SleepstageW = [];
    annotations.SleepstageN1 = [];
    annotations.Sleepstage1 = [];
    annotations.SleepstageN2 = [];
    annotations.Sleepstage2 = [];
    annotations.SleepstageN3 = [];
    annotations.Sleepstage3 = [];
    annotations.Sleepstage4 = [];
    annotations.SleepstageR = [];
    for i = 1:length(edfAnnotations.label)
        try 
            name = edfAnnotations.label{i};
            aux = strfind(edfAnnotations.label{i}, '@@');
            if ~isempty(aux) || strncmp(name, 'Sleep stage', 11)
                if strncmp(name, 'Sleep stage', 11)
                    aux = length(name) + 1;
                end
                name = regexprep(name(1:aux - 1), '[ (),_?]', '');
                if isfield(annotations, name)
                    last = length(annotations.(name));
                else
                    last = 0;
                end
                start = edfAnnotations.offset(i) - deltaTime;
                duration = edfAnnotations.duration(i);
                if strncmp(name, 'Sleepstage', 10) 
                    if edfAnnotations.duration(i) < 10
                        logger.warning('Skipping annotation: %s (%d +%d s)', name, start, duration);
                        continue
                    end
                    start = round(start / 30) * 30;
                    duration = round(duration / 30) * 30;
                end
                annotations.(name)(last + 1).start = start;
                annotations.(name)(last + 1).duration = duration;
            end 
        catch
            logger.warning('Strange annotation: %s', name);
        end
    end
end

function time = stringToSeconds(hourString)
    hours = str2double(hourString(1:2));
    minutes = str2double(hourString(4:5));
    seconds = str2double(hourString(7:8));
    
    time = hours * 60 * 60 + minutes * 60 + seconds;
end