%FINDALPHAAROUSALS Find alpha arousals using the initial events
%   AR = FINDALPHAAROUSALS(EEG, EMG) returns the alpha arousals detected in
%   the EEG signal studying the initial EEG.EVENTS and using the EMG signals
%   to improve the detection. 
%
%   [AR, F] = FINDALPHAAROUSALS(EEG, EMG) also returns the list of events
%   that are discarded as arousal with the reason in F.INFO
%
%   [AR, F, S] = FINDALPHAAROUSALS(EEG, EMG) returns an empty list. 

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

function [arousals, removed, spindles] = findAlphaArousals(eeg, emg)
    logger = logging.Logger.getLogger('FindAlphaArousals');
    
    logger.info('Analysis of alpha events...');

    alphaEventsArray = eeg.events.alpha;
    
    [eventsStart, eventsEnd] = getEventIndexes(alphaEventsArray);
    
    % Malloc structure
    events(length(eventsStart)).startSample = 0;
    
    for i = 1:length(eventsStart)
        events(i).startIdx = eventsStart(i);
        events(i).endIdx = eventsEnd(i);
        events(i).start = (eventsStart(i) - 1) / eeg.eventsRate;
        events(i).duration = (eventsEnd(i) - eventsStart(i)) / eeg.eventsRate;
        events(i).startSample = round(events(i).start * eeg.rate);
        events(i).endSample = events(i).startSample + round(events(i).duration * eeg.rate);
        events(i).info = {['Alpha event: ' num2str(events(i).start) ' + ' num2str(events(i).duration) ' s']};
    end
    
    % Remove events happening at the beginning and the end to avoid future errors
    events([events.start] < 60) = [];
    events([events.startSample] > length(eeg.raw) - 60 * eeg.rate) = [];
    % Remove events too short
    events([events.startIdx] == [events.endIdx]) = [];

    global updateStats;
    
    updateStats.betaEventsUpdate = 0;
    updateStats.updateWithOverlappingEmg = 0;
    updateStats.ampUpdate = 0;
    updateStats.updateWithNonContinuousOverlappingEmg = 0;
    updateStats.updateWithNonOverlappingEmg = 0;
    updateStats.powerUpdate = 0;
    
    spindleEvents = [];
    arousalEvents = false(1, length(events));
    for i = 1:length(events)
        event = events(i);
        
        originalDuration = event.duration;
        
        % We consider updating the event only if a certain minimum original duration
        % is reached
        if (originalDuration > 1)
            logger.finer('Studing event from %f (+ %f) s', event.start, event.duration);
            event = updateEvent(event, eeg, emg, eeg.energies.alpha);
        
            emgActivity = checkEmgActivity(event, emg);
            if ~emgActivity
                logger.finer('Event without EMG activity');
                event.info{end + 1} = 'No EMG activity';
            else
                logger.finer('Event with EMG activity');
                event.info{end + 1} = 'EMG activity';
            end
        end
        
        validArousal = 1;
        
        if event.duration < 4
            validArousal = 0; 
            event.info{end + 1} = 'Removed: too short';
        end
        
        timeLimit = 20;
        if event.duration > timeLimit
            endSecond = event.start + event.duration;
            epoch1 = ceil(event.start / 30); 
            epoch2 = ceil(endSecond / 30);
            if epoch1 ~= epoch2
                firstEpoch = epoch1 * 30 - event.start;
                secondEpoch = endSecond - epoch1 * 30;
                if firstEpoch > timeLimit || secondEpoch > timeLimit
                    validArousal = 0; 
                    event.info{end + 1} = 'Removed: too long';
                end
            else
                validArousal = 0;
                event.info{end + 1} = 'Removed: too long';
            end
        end
                
        if validArousal && ~emgActivity
            validArousal = 0;
            event.info{end + 1} = 'Removed: no EMG';
            logger.finest('Alpha removed: no EMG (%f + %f)', event.start, event.duration);
        end
        
        arousalEvents(i) = validArousal;
        
        events(i) = event;
    end
    
    logger.info('BetaEvents updates: %d', updateStats.betaEventsUpdate);
    logger.info('Power updates: %d', updateStats.powerUpdate);
    logger.info('Amplitude updates: %d', updateStats.ampUpdate);
    logger.info('Overlapping EMG updates: %d', updateStats.updateWithOverlappingEmg);
    logger.info('Non-continuous overlapping EMG updates: %d', updateStats.updateWithNonContinuousOverlappingEmg);
    logger.info('Non overlapping EMG updates: %d', updateStats.updateWithNonOverlappingEmg);
    
    arousals = events(arousalEvents);
    removed = events(~arousalEvents);
    spindles = spindleEvents;
end

function emgActivity = checkEmgActivity(event, emg)
    minEMGDuration = 1; % In seconds

    eventPosition = round(event.start * emg.eventsRate):...
            round((event.start + event.duration) * emg.eventsRate);
        
    emgActivity = sum(emg.events(eventPosition)) > (minEMGDuration * emg.eventsRate);
end

function event = updateEvent(event, eeg, emg, power)
    global updateStats;
            
    oldDuration = event.duration;
    event = updateWithPowerAlpha(event, eeg, power);
    if event.duration > oldDuration, updateStats.powerUpdate = updateStats.powerUpdate + 1; end;

    oldDuration = event.duration;
    event = updateWithOverlappingEmg(event, emg);
    if event.duration > oldDuration, updateStats.updateWithOverlappingEmg = updateStats.updateWithOverlappingEmg + 1; end;

    oldDuration = event.duration;
    event = updateWithNonContinuousOverlappingEmg(event, emg);
    if event.duration > oldDuration, updateStats.updateWithNonContinuousOverlappingEmg = updateStats.updateWithNonContinuousOverlappingEmg + 1; end;

    oldDuration = event.duration;
    event = updateWithNonOverlappingEmg(event, emg);
    if event.duration > oldDuration, updateStats.updateWithNonOverlappingEmg = updateStats.updateWithNonOverlappingEmg + 1; end;

    oldDuration = event.duration;
    event = updateWithAmplitude(event, eeg);
    if event.duration > oldDuration, updateStats.ampUpdate = updateStats.ampUpdate + 1; end;
        
end
