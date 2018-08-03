%updateWithNonOverlappingEmg Update the event using the emg
%signal. 
%   EVENT = updateWithNonOverlappingEmg(EVENT, EEG) returns the 
%   event with the info and the duration updated

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

function event = updateWithNonOverlappingEmg(event, emg)    
    logger = logging.Logger.getLogger('updateWithNonOverlappingEmg');

    limit = Configuration.instance().updateWithFutureEmgLimit;
    windowSize = selectWindowSize(emg.rate);
    
    threshold = Configuration.instance().updateWithFutureEmgThreshold;

    amplitudeWindowSize = selectAmplitudeWindowSize(emg.rate);
    amplitudeWindowSamples = round(amplitudeWindowSize * emg.rate);
    
    ampWindowsWindow = floor(windowSize / amplitudeWindowSize);
      
    prevWindow = abs(emg.raw(event.startSample - 10 * emg.rate:event.startSample - 1 * emg.rate - 1));
    prevWindow = prevWindow(1:floor(length(prevWindow) / amplitudeWindowSamples) * amplitudeWindowSamples);
    % Median amplitude for amplitudeWindowSize seconds window
    previousAmp = median(max(reshape(prevWindow, [amplitudeWindowSamples, length(prevWindow) / amplitudeWindowSamples])));
    previousAmp = max(previousAmp, 15);
    
    % Peak-to-peak amplitude for the next 'limit' seconds in amplitudeWindowSize s windows
    if (event.endSample + limit * emg.rate - 1 > length(emg.raw))
        return;
    end
    futureWindow = abs(emg.raw(event.endSample:event.endSample + limit * emg.rate - 1));
    futureWindow = futureWindow(1:floor(length(futureWindow) / amplitudeWindowSamples) * amplitudeWindowSamples);
    futureAmp = max(reshape(futureWindow, [amplitudeWindowSamples, length(futureWindow) / amplitudeWindowSamples]));
    
    numWindows = floor(length(futureAmp) / ampWindowsWindow);
    
    bigEmg = 0;
    smallEmg = 0;
    newEnd = 0;
    
    for i = 1:numWindows
        amplitude = median(futureAmp((i - 1) * ampWindowsWindow + 1:i*ampWindowsWindow));
        
        if amplitude > previousAmp * threshold
            bigEmg = bigEmg + 1;
            smallEmg = 0;
        else
            smallEmg = smallEmg + 1;
            if bigEmg >= 2
                newEnd = i - 1;
            end
            if i - 1 - newEnd >= 4 / windowSize
                break
            end
            bigEmg = 0;
        end
    end
    
    if i == numWindows && bigEmg > 0
        newEnd = numWindows;
    end
    
    if newEnd > 0
        duration = event.duration + newEnd * windowSize;

        if duration > 15 
            logger.info('long event starting in %f + %f s (previous: %f s)', event.start, duration, event.duration);
        end
                
        event.info{end + 1} = ['EMG future + ' num2str(duration - event.duration) ' s'];
        event.endSample = event.startSample + round(duration * emg.rate);
        event.duration = duration;  
    end
end

function time = selectWindowSize(rate) 
    if rate == 125 || rate == 250
        time = 0.48;
    else
        time = 0.5;
    end
end

function time = selectAmplitudeWindowSize(rate) 
    if rate == 125 || rate == 250
        time = 0.12;
    else
        time = 0.125;
    end
end