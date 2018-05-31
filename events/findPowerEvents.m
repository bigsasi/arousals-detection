%FINDPOWEREVENTS Find power based events in the signal
%   signal = findPowerEvents(signal) returns the original signal adding
%   the following fields to the structure:
%      signal.energies: list of vectors with the power used for finding the
%      events in each power band.
%      signal.baseline: list of vectors with the baseline used for finding
%      the events in each power band.
%      signal.events: the events founded. 
%      signal.eventsRate: the events rate.

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

function signal = findPowerEvents(signal)
    windowStep = Configuration.instance().powerEventsWindowStep;
    windowSize = Configuration.instance().powerEventsWindowSize;
    numWindows = round(0.5 * windowSize / windowStep);
    backgroundTime = Configuration.instance().powerEventsBackgroundTime;

    [alpha, emg] = signalPower(signal.raw, windowSize, signal.rate, windowStep, [8 16], [12 signal.rate / 2]);

    power.alpha = alpha;
    power.emg = emg;

    power.alpha(numWindows:end) = power.alpha(1:end - numWindows + 1);
    power.emg(numWindows:end) = power.emg(1:end - numWindows + 1);

    baselines.alpha = baseline(backgroundTime, windowStep, power.alpha);
    baselines.emg = baseline(backgroundTime, windowStep, power.emg);
        
    events.alpha = getSignalEvents(power.alpha, Configuration.instance().powerEventsAlphaThreshold, 0, windowStep, baselines.alpha);
    events.emg = getSignalEvents(power.emg, Configuration.instance().powerEventsBetaThreshold, 0, windowStep, baselines.emg);
    
    if Configuration.instance().updateWithAlphaEvents
        [events.emg, updates] = updateAlphaBetaOverlapping(events.emg, events.alpha, Configuration.instance().updateWithAlphaEventsUnionTime / windowStep);
        while (updates > 0)
            [events.emg, updates] = updateAlphaBetaOverlapping(events.emg, events.alpha, Configuration.instance().updateWithAlphaEventsUnionTime / windowStep);
        end
    end
    if Configuration.instance().updateWithBetaEvents
        [events.alpha, updates] = updateAlphaBetaOverlapping(events.alpha, events.emg, Configuration.instance().updateWithBetaEventsUnionTime / windowStep);
        while (updates > 0)
            [events.alpha, updates] = updateAlphaBetaOverlapping(events.alpha, events.emg, Configuration.instance().updateWithBetaEventsUnionTime / windowStep);
        end
    end
    
    % Merge events separated by less than a certain distance
    events.alpha = mergeEventsCloserThan(events.alpha, Configuration.instance().powerEventsBetaUnion, windowStep);
    events.emg = mergeEventsCloserThan(events.emg, Configuration.instance().powerEventsAlphaUnion, windowStep);

    signal.events = events;
    signal.eventsRate = 1 / windowStep;
    signal.energies = power;
    signal.baselines = baselines;
end
