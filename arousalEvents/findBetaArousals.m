%FINDBETAAROUSALS Find beta arousals using the initial events
%   AR = FINDBETAAROUSALS(EEG, EMG) returns the beta arousals detected in
%   the EEG signal studying the initial EEG.EVENTS and using the EMG signals
%   to improve the detection. 
%
%   [AR, F] = FINDBETAAROUSALS(EEG, EMG) also returns the list of events
%   that are discarded as arousal with the reason in F.INFO
%
%   [AR, F, S] = FINDBETAAROUSALS(EEG, EMG) also returns the list of
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

function [arousals, removed, spindles] = findBetaArousals(eeg, emg)
    global updateStats;
    logger = logging.Logger.getLogger('FindBetaArousals');
    
    betaEventsArray = eeg.events.emg;
    
    [eventsStart, eventsEnd] = getEventIndexes(betaEventsArray);
    
    % Malloc structure
    events(length(eventsStart)).startSample = 0;
    
    for i = 1:length(eventsStart)
        events(i).startIdx = eventsStart(i);
        events(i).endIdx = eventsEnd(i);
        events(i).start = (eventsStart(i) - 1) / eeg.eventsRate;
        events(i).duration = (eventsEnd(i) - eventsStart(i)) / eeg.eventsRate;
        events(i).startSample = round(events(i).start * eeg.rate);
        events(i).endSample = round((events(i).start + events(i).duration) * eeg.rate);
        events(i).info = {['Beta event: ' num2str(events(i).start) ' + ' num2str(events(i).duration) ' s']};
    end
    
    % Remove events happening at the beginning to avoid future errors
    events([events.start] < 60) = [];
    events([events.startSample] > length(eeg.raw) - 60 * eeg.rate) = [];
    % Remove events too short
    events([events.startIdx] == [events.endIdx]) = [];

    logger.info('Starting with %d events', length(events)); 
       
    eegSmooth = eeg;
    eegSmooth.raw = polymanBandPassFilter(eeg.raw, eeg.rate, 1, 1, 32);
    cycles = signalCycles(eegSmooth.raw);
    
    updateStats.alphaEventsUpdate = 0;
    updateStats.updateWithOverlappingEmg = 0;
    updateStats.ampUpdate = 0;
    updateStats.updateWithNonContinuousOverlappingEmg = 0;
    updateStats.updateWithNonOverlappingEmg = 0;
    updateStats.powerUpdate = 0;
    
    spindleEvents = [];
    arousalEvents = false(1, length(events));
    for i = 1:length(events)
        event = events(i);
        
        logger.finer('Studing event from %f (+ %f) s', event.start, event.duration);
        event = updateEvent(event, eeg, emg, eeg.energies.emg);
                
        emgActivity = checkEmgActivity(event, emg);
        if ~emgActivity
            logger.finer('Event without EMG activity');
            event.info{end + 1} = 'No EMG activity';
        else
            logger.finer('Event with EMG activity');
            event.info{end + 1} = 'EMG activity';
        end
        event.endIdx = event.startIdx + floor(event.duration * eeg.eventsRate);

        validArousal = 1;
        
        if event.duration < 3, 
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
        
        % Remove spindles if there is no EMG activity
        if validArousal && ~emgActivity
            [event, spindles, spindleFree] = removeSpindles(event, eegSmooth, cycles);
            
            if ~spindleFree, 
                removedSpindlesString =  ['Removed spindles in ' num2str(event.start) ' + ' num2str(event.duration) ': '];
                for j = 1:length(spindles)
                    removedSpindlesString = [removedSpindlesString num2str(spindles(j).start) ' (' num2str(spindles(j).duration) ') '];
                end
                logger.info(removedSpindlesString);
                validArousal = 0;
                event.info{end + 1} = 'Removed: spindles';
            end;
            if ~isempty(spindles), spindleEvents = [spindleEvents spindles]; end;
        end
        
        arousalEvents(i) = validArousal;
        
        events(i) = event;
    end
    
    logger.info('AlphaEvents updates: %d', updateStats.alphaEventsUpdate);
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
    eventPosition = round(event.start * emg.eventsRate):...
            round((event.start + event.duration) * emg.eventsRate);
    
    emgActivity = any(emg.events(eventPosition));        
end


function event = updateEvent(event, eeg, emg, power)
    global updateStats;
    
    oldDuration = event.duration;
    event = updateWithPowerBeta(event, eeg, power);
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
