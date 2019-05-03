%UPDATEWITHBETAEVENTS Check if the event overlaps with beta events on a 
%certain neighborhood
%   EVENT = updateWithBetaEvents(EVENT, EEG) returns the updated event if it
%   overlaps with another beta events

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

function event = updateWithBetaEvents(event, eeg)
    unionTime = round(Configuration.updateWithBetaEventsUnionTime * eeg.eventsRate);
    
    searchInterval = eeg.events.emg(event.startIdx-unionTime:event.endIdx+unionTime);
    
    if any(searchInterval)
        oldDuration = event.duration;
        
        iniIdx = find(searchInterval,1,'first');
        endIdx = find(searchInterval,1,'last');
        
        % Check the real start of the first overlapping alpha event
        if (iniIdx == 1)
            k = 0;
            while eeg.events.emg(event.startIdx-unionTime-k-1)
                k = k + 1;
            end
            event.startIdx = event.startIdx-unionTime-k;
        end
        
        % Check the real end of the last overlapping alpha event
        if (endIdx == length(searchInterval))
            k = 0;
            while eeg.events.emg(event.endIdx+unionTime+k+1)
                k = k + 1;
            end
            event.endIdx = event.endIdx+unionTime+k;
        end    
        
        % Recompute final fields for the resulting event
        if ne(oldDuration, (event.endIdx - event.startIdx) / eeg.eventsRate) 
            event.start = (event.startIdx - 1) / eeg.eventsRate;
            event.duration = (event.endIdx - event.startIdx) / eeg.eventsRate;
            event.startSample = round(event.start * eeg.rate);
            event.endSample = round((event.start + event.duration) * eeg.rate);
            event.info{end + 1} = ['BetaEvent update +' num2str(event.duration - oldDuration) ' s'];
        end
    end
end

