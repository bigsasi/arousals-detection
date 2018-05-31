%REMOVEAWAKEEVENTS Remove arousals during W stage using an hypnogram
%   EVENTS = REMOVEAWAKEEVENTS(EVENTS, HYPNOGRAM) returns the events
%   happening during an epoch not scored as W
% 
%   [EVENTS, DISCARDED] = REMOVEAWAKEEVENTS(EVENTS, HYPNOGRAM) also returns
%   the list of events happening during W epochs.

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

function [validEvents, falseEvents] = removeAwakeEvents(events, hypnogram)
    awaken = false(1, length(events));
    
    for i = 1:length(events)
        from = events(i).start;
        to = from + events(i).duration;
        
        fromEpoch = min(length(hypnogram), floor(from / 30) + 1);
        toEpoch = min(length(hypnogram), floor(to / 30) + 1);
        
        awakeFrom = hypnogram(fromEpoch) == 0;
        awakeTo = hypnogram(toEpoch) == 0;
        
        awaken(i) = awakeFrom || awakeTo;
    end
    
    validEvents = events(~awaken);
    falseEvents = events(awaken);
end

