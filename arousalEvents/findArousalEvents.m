%FINDAROUSALEVENTS Find arousals using the initial events
%   AR = FINDAROUSALEVENTS(EEG, EMG, HYP) returns the arousals detected in
%   the EEG signal studying the initial EEG.EVENTS and using the EMG signals
%   and the hypnogram HYP to improve the detection.
%
%   [AR, F] = FINDAROUSALEVENTS(EEG, ...) also returns the list of events
%   that are discarded as arousal with the reason in F.INFO
%
%   [AR, F, S] = FINDAROUSALEVENTS(EEG, ...) also returns the list of
%   spindles found in the initial events. 

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

function [arousals, removedArousals, spindles] = findArousalEvents(eeg, emg, hypnogram)
    logger = logging.Logger.getLogger('FindArousalEvent');
    
    betaArousals = [];
    removedBeta = [];
    betaSpindles = [];
    if any(eeg.events.emg)
        [betaArousals, removedBeta, betaSpindles] = findBetaArousals(eeg, emg);
        for i = 1:length(betaArousals)
            betaArousals(i).source = 'beta';
        end
        arousals = betaArousals;
    
        eeg = removeBetaArousals(eeg, betaArousals);
    end
    
    alphaSpindles = [];
    removedAlpha = [];
    
    if any(eeg.events.alpha)
        [alphaArousals, removedAlpha, alphaSpindles] = findAlphaArousals(eeg, emg);
        for i = 1:length(alphaArousals)
            alphaArousals(i).source = 'alpha';
        end

        logger.info('Detected: %d alpha arousals', length(alphaArousals));

        arousals = mergeEvents(alphaArousals, betaArousals);
    end
        
    removedArousals = [];
    if (~isempty(arousals))
        [arousals, removedAwake] = removeAwakeEvents(arousals, hypnogram);
        for i = 1:length(removedAwake)
            removedAwake(i).info{end + 1} = 'Removed: awake';
        end
        
        logger.info('Removed %d arousals for being awake', length(removedAwake));
        
        removedArousals = removedAwake;
    end
            
    removedRem = [];
    if (~isempty(arousals))
%         [arousals, removedRem] = removeRemEvents(arousals, emg);
        [arousals, removedRem] = removeCheckREMeventsByHypnogram(arousals, emg, hypnogram);
    
        logger.info('Removed %d arousals in REM', length(removedRem));
    end
    
    if (~isempty(arousals))
        arousals = removeIfClose(arousals, 10);
    end    
    
    arousals = removeIncluded(arousals);
    
    removedArousals = mergeEvents(removedArousals, removedRem);
    removedArousals = mergeEvents(removedArousals, removedBeta);
    removedArousals = mergeEvents(removedArousals, removedAlpha);

    spindles = mergeEvents(betaSpindles, alphaSpindles);
end

function events = removeIncluded(events)
    remove = false(length(events), 1);
    for i = 1:length(events) - 1
        if events(i).start + events(i).duration > events(i + 1).start
            remove(i + 1) = 1;
        end
    end
    events(remove) = [];
end

function [validEvents, falseEvents] = removeRemEvents(events, emg)
    
    windowSize = 0.1;
    numWindows = 30 / windowSize;
    oneSecEmg = signalAmplitude(emg.raw, emg.rate, windowSize, windowSize);
    lengthWindows = round(length(oneSecEmg) / numWindows);
    if lengthWindows * numWindows < length(oneSecEmg)
        oneSecEmg(lengthWindows * numWindows + 1:end) = [];
    elseif lengthWindows * numWindows > length(oneSecEmg)
        oneSecEmg(end + 1:lengthWindows * numWindows) = oneSecEmg(end);
    end
    epochAmplitudes = median(reshape(oneSecEmg, [numWindows, lengthWindows]));
    tmp = sort(epochAmplitudes);
    remLimit = tmp(round(0.2 * length(epochAmplitudes)));

    remove = false(1, length(events));
    
    for i = 1:length(events)
        from = events(i).start;
        to = from + events(i).duration;
        
        fromEpoch = min(length(epochAmplitudes), floor(from / 30) + 1);
        toEpoch = min(length(epochAmplitudes), floor(to / 30) + 1);
        
        isRem = any(epochAmplitudes([fromEpoch toEpoch]) < remLimit);
        
        remove(i) = isRem && ~checkEmgActivity(events(i), emg);
    end

    validEvents = events(~remove);
    falseEvents = events(remove);
    for i = 1:length(falseEvents)
        falseEvents(i).info{end + 1} = 'Removed: rem';
    end
end

function [validEvents, falseEvents] = removeCheckREMeventsByHypnogram(events, emg, hypnogram)
    
    remove = false(1, length(events));
    
    for i = 1:length(events)
        from = events(i).start;
        to = from + events(i).duration;
                       
        fromEpoch = min(length(hypnogram), floor(from / 30) + 1);
        toEpoch = min(length(hypnogram), floor(to / 30) + 1);
                        
        remFrom = hypnogram(fromEpoch) == 5;
        remTo = hypnogram(toEpoch) == 5;
        
        isREM = remFrom && remTo;
        
        remove(i) = isREM && ~checkEmgActivity(events(i), emg);
    end

    validEvents = events(~remove);
    falseEvents = events(remove);
    for i = 1:length(falseEvents)
        falseEvents(i).info{end + 1} = 'Removed: rem';
    end
end

function emgActivity = checkEmgActivity(event, emg)
    eventPosition = round(event.start * emg.eventsRate):...
            round((event.start + event.duration) * emg.eventsRate);
    
    emgActivity = any(emg.events(eventPosition));        
end


function eeg = removeBetaArousals(eeg, arousals)
    for i = 1:length(arousals)
        arousalStart = arousals(i).start * eeg.eventsRate;
        arousalEnd = (arousals(i).start + arousals(i).duration) * eeg.eventsRate;
        
        arousalStart = round(arousalStart) - 10 * eeg.eventsRate;
        arousalEnd = round(arousalEnd);
        
        eeg.events.alpha(arousalStart:arousalEnd) = 0;
    end
end

function validEvents = removeIfClose(events, time)

    keep = true(length(events), 1);
    for i = 2:length(events)
        previousEnd = events(i - 1).start + events(i - 1).duration;
        if keep(i - 1) && events(i).start <= previousEnd + time
            keep(i) = 0;
        end
    end
    validEvents = events(keep);
end

