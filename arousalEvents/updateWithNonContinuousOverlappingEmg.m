%updateWithNonContinuousOverlappingEmg Update the event using the emg
%signal. 
%   EVENT = updateWithNonContinuousOverlappingEmg(EVENT, EEG) returns the 
%   event with the info and the duration updated

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

function event = updateWithNonContinuousOverlappingEmg(event, emg)
    logger = logging.Logger.getLogger('updateWithNonContinuousOverlappingEmg');

    cutSize = selectCutSize(emg.rate);

    previousWindowAxe = event.startSample - 3 * emg.rate:event.startSample - 1;
    currentWindowAxe = previousWindowAxe;
    
    limit = floor(60 / cutSize);
    eventLimit = floor(event.duration / cutSize);
    
    if event.duration < cutSize
        return
    end

    % Mean difference between two rectified emg windows
    previousDiff = mean(abs(emg.raw(previousWindowAxe - cutSize * emg.rate) - emg.raw(previousWindowAxe)));
    
    eventDiff = zeros(1, eventLimit);
    for idx = 1:eventLimit
        currentWindowAxe = currentWindowAxe + cutSize * emg.rate;
        eventDiff(idx) = mean(abs(emg.raw(currentWindowAxe) - emg.raw(previousWindowAxe)));
    end
    
    % Values to decide if the update is necessary
    eventDiffMean = mean(eventDiff);
    eventDiffStd = std(eventDiff);
    
    if all(eventDiff < 2 * previousDiff)
        return
    end
    
    postDiff = nan(1, limit);
    for idx = eventLimit + 1:limit
        currentWindowAxe = currentWindowAxe + cutSize * emg.rate;
        postDiff(idx) = mean(abs(emg.raw(currentWindowAxe) - emg.raw(previousWindowAxe)));
        
        % Update the event while the conditions hold
        if postDiff(idx) < eventDiffMean - eventDiffStd
            break
        end
        if sum(postDiff < 1.1 * (eventDiffMean - eventDiffStd)) > 2
            break
        end
    end

    duration = event.duration + (idx - eventLimit - 6) * cutSize;
        
    if duration > event.duration
        if duration > 15
            logger.info('discarding update in %f + %f s (previous: %f s)', event.start, duration, event.duration);
            return
        end
        event.info{end + 1} = ['EMG similarity +' num2str(duration - event.duration) ' s'];
        event.endSample = event.startSample + round(duration * emg.rate);
        event.duration = duration;    
    end 
end


function time = selectCutSize(rate) 
    if rate == 125 || rate == 250
        time = 0.52;
    else
        time = 0.5;
    end
end