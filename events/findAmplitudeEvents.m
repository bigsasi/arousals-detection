%FINDAMPLITUDEEVENTS Find amplitude based events in the signal
%   signal = findAmplitudeEvents(signal) returns the original signal adding
%   the following fields to the structure:
%      signal.amplitude: the amplitude values used to obtain the events.
%      signal.baselineAmplitude: the amplitude baseline used to obtain the
%      events.
%      signal.meanAmplitude: the mean amplitude of the signal.
%      signal.events: the events founded. 
%      signal.eventsRate: the events rate.

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

function signal = findAmplitudeEvents(signal)
    threshold = Configuration.instance().amplitudeEventsThreshold;
    windowStep = Configuration.instance().amplitudeEventsWindowStep;
    windowSize = Configuration.instance().amplitudeEventsWindowSize;
    backgroundTime = Configuration.instance().amplitudeEventsBackgroundTime;
    
    numWindows = round(0.5 * windowSize / windowStep);

    amplitude = signalAmplitude(signal.raw, signal.rate, windowSize, windowStep);
    baselineAmplitude = baseline(backgroundTime, windowStep, amplitude);
    meanAmplitude = ones(1, length(amplitude)) * mean(amplitude);

    amplitude(numWindows:end) = amplitude(1:end - numWindows + 1);
    baselineAmplitude(numWindows:end) = baselineAmplitude(1:end - numWindows + 1);

    events = getSignalEvents(amplitude, threshold, Configuration.instance().amplitudeEventsMergeNear, windowStep, baselineAmplitude);
    
    % Clean really small events
    for i = 1:length(events) - 3
        if all(events(i:i+2) == [0; 1; 0])
            events(i + 1) = 0;
        end
    end
            
    signal.amplitude = amplitude;
    signal.baselineAmplitude = baselineAmplitude;
    signal.meanAmplitude = meanAmplitude;
    signal.events = events;
    signal.eventsRate = 1 / windowStep;
end