%UPDATEWITHPOWERALPHA Update the event using alpha power values
%   EVENT = updateWithPowerAlpha(EVENT, EEG, POWER) returns the  event with 
%   the info and the duration updated

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

function event = updateWithPowerAlpha(event, eeg, power)
    logger = logging.Logger.getLogger('UpdateWithPowerAlpha');

    powerRate = eeg.eventsRate;
    threshold = Configuration.instance().updateWithPowerAlphaThreshold;
    
    limitTime = Configuration.instance().updateWithPowerAlphaLimitTime;
    increaseTime = Configuration.instance().updateWithPowerAlphaIncreaseTime;
        
    previousSeconds = Configuration.instance().updateWithPowerAlphaPreviousSeconds;
       
    previousPowerSamples = previousSeconds * powerRate;
    
    endIdx = event.startIdx + round(event.duration * powerRate);

    previousWindow = event.startIdx - previousPowerSamples:event.startIdx - 1;
    
    previousPower = mean(power(previousWindow));
    
    % Discard if it is just a strange peak
    previousSignalSamples = 10 * eeg.rate;
    previousWindow = event.startSample - previousSignalSamples:event.startSample - 1;
    eventWindow = event.startSample:event.endSample;
    previousAmp = max(abs(eeg.raw(previousWindow)));
    eventAmp = max(abs(eeg.raw(eventWindow)));
    if eventAmp > 8 * previousAmp
        return
    end

    
    lastWindow = endIdx;
    newEnd = endIdx;
    while newEnd < endIdx + 30 * powerRate
        nextPower = mean(power(lastWindow + 1:lastWindow + increaseTime * powerRate));
        if nextPower > threshold * previousPower
            newEnd = lastWindow + increaseTime * powerRate;
        end
                
        lastWindow = lastWindow + increaseTime * powerRate;
        
        if lastWindow >= newEnd + limitTime * powerRate
            break
        end
    end
          
    if newEnd > endIdx
        duration = (newEnd - endIdx) / powerRate;
        if duration > 15
            logger.info('discarding update in %f + %f s (previous: %f s)', event.start, duration, event.duration);
            return
        end
        event.info{end + 1} = ['Power update + ' num2str(duration) ' s'];
        event.duration = event.duration + duration;
        event.endSample = event.startSample + round(duration * eeg.rate);
    end       
end