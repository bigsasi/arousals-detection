%LOADSHHSANNOTATIONS load annotations from a xml file. xml structure should follow 
% shhs xml files. This is, events should be annotated following the structure
% CMPStudyConfig.ScoredEvents.ScoredEvent
%   filename    complete file name (with path) 
%   annotations struct with a field for each event type. Field names are the 
%               ones used in the xml files, without white spaces.

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

function [annotations, hypnogram] = loadXmlAnnotations(file)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    annotations = [];
    xml = xml_parseany(file);
    totalEvents = length(xml.ScoredEvents.ScoredEvent);
    for i = 1:totalEvents
        currentEvent = xml.ScoredEvents.ScoredEvent(i);
        name = regexprep(currentEvent.Name, '[ ()]', '');
        if isfield(annotations, name)
            last = length(annotations.(name));
        else
            last = 0;
        end
        annotations.(name)(last + 1).start = currentEvent.Start;
        annotations.(name)(last + 1).duration = currentEvent.Duration;
    end
    
    hypnogram = zeros(length(xml.SleepStages.SleepStage), 1);
    for i = 1:length(xml.SleepStages.SleepStage)
        hypnogram(i) = xml.SleepStages.SleepStage{i};
    end
end

