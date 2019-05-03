%REMOVESPINDLES Detect spindles in an event
%   EVENT = REMOVESPINDLES(EVENT, EEG, DIFFS) returns the EVENT adding
%   info about the founded spindles, if any. DIFFS should be the aproximate
%   second derivative of the EEG signal. 
%
%   [EVENT, SP] = REMOVESPINDLES(EVENT, EEG, DIFFS) returns the list of
%   spindles founded inside the event.
%
%   [EVENT, SP, SAVE] = REMOVESPINDLES(EVENT, EEG, DIFFS) returns a boolean
%   true if the event is still valid or false if it should be discarded
%   because of the presence of spindles.

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

function [event, spindles, saveEvent] = removeSpindles(event, smooth, cycles)
    spindles = [];
    saveEvent = 1;
    
    % Find spindles inside the event
    [foundSpindles] = findSpindle(event, smooth, cycles);

    if ~isempty(foundSpindles)
        % Remove the spindles from the event
        [remainingArousal] = remainingEvent(event, foundSpindles);       

        spindles = foundSpindles([foundSpindles.duration] > 0);

        % Save the spindles if a remaining portion of the event lasts more
        % than 3 seconds
        idx = find([remainingArousal.duration] >= 3, 1);

        if isempty(idx)
            saveEvent = 0;
        else
            for j = 1:length(foundSpindles)
                event.info{end + 1} = ['Spindle: ' num2str(foundSpindles(j).start) ' + ' num2str(foundSpindles(j).duration)];
            end
        end
    end
end

function events = remainingEvent(event, toRemove)
    events(1:length(toRemove) + 1) = event;
    
    for i = 1:length(toRemove) 
        events(i).duration = toRemove(i).start - events(i).start;
        events(i + 1).start = toRemove(i).start + toRemove(i).duration;
    end
    
    events(end).duration = event.duration + event.start - events(end).start;
end

function numCycles = countCycles(signal)
      numCycles = sum(signal > 0);
end

function [spindles] = findSpindle(event, signal, cycles)
    
    spindles(1) = event;
    spindles(1).duration = -1;
    lastSpindle = 1;
    
    firstSample = round((event.start - 0.5) * signal.rate);
    lastSample = round((event.duration + event.start) * signal.rate);
    
    nextSpindle = round(signal.rate * 0.125);
    deltaTime = round(signal.rate * 0.125);
    windowSize = round(signal.rate * 0.5);
    
    i = firstSample;
    while i < lastSample - windowSize
        % -2 because cycles is the second derivative
        numCycles = countCycles(cycles(i:i + windowSize - 2));
        
        freq = numCycles / 0.5;
        j = i + windowSize;
        
        numCycles = countCycles(cycles(i:j - 2));
        while (freq >= 12) && (freq <= 14) && (j + deltaTime < length(signal.raw))
            j = j + deltaTime;
%             numCycles = countCycles(cycles(i:j - 2));
            numCycles = numCycles + sum(cycles(j - deltaTime - 1: j - 2) > 0);
            freq = numCycles / ((j - i) / signal.rate);
        end
        
        isSpindle = (j - deltaTime - i) > 0.5 * signal.rate;
        
        if isSpindle
            spindles(lastSpindle) = event;
            spindles(lastSpindle).start = i / signal.rate;
            spindles(lastSpindle).duration = (j - deltaTime - i) / signal.rate;
            lastSpindle = lastSpindle + 1;
            i = j + nextSpindle;
        else
            i = i + deltaTime;
        end
        
    end
    
    if spindles(1).duration == -1
        spindles = [];
    end
end