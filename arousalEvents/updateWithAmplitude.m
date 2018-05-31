%UPDATEWITHAMPLITUDE update the event using the voltage values from the eeg
%signal
%   EVENT = updateWithAmplitude(EVENT, EEG) returns the event with the info
%   and the duration updated

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

function event = updateWithAmplitude(event, eeg)
    threshold = Configuration.instance().updateWithAmplitudeThreshold;
        
    startSample = event.startSample;
    endSample = min(length(eeg.raw), event.endSample);
    
    prevSample = startSample - eeg.eventsRate * eeg.rate;
    allAmps = signalAmplitude(eeg.raw(prevSample:startSample), eeg.rate, 1, 1);
    
    % mean 1 second window peak-to-peak amplitude for the previous 5 seconds 
    prevAmp = mean(allAmps(allAmps > 0));

    % Event peak-to-peak amplitude
    eventAmp = max(eeg.raw(startSample:endSample)) - min(eeg.raw(startSample:endSample));

    if eventAmp > threshold * prevAmp
        for sample = endSample + eeg.rate:eeg.rate:length(eeg.raw) - eeg.rate
            % Add 1 second window to the event
            allAmps = signalAmplitude(eeg.raw(sample:sample + eeg.rate), eeg.rate, 1, 1);
            % Peak-to-peak amplitude of the added window
            eventAmp = mean(allAmps(allAmps > 0));

            if eventAmp > threshold * prevAmp
                endSample = sample;
            else
                % Finish the update
                break
            end
        end
    end

    duration = (endSample - startSample) / eeg.rate;
    if duration > event.duration + 1
        event.info{end + 1} = ['Amplitude update +' num2str(duration - event.duration) ' s'];
        event.duration = duration;
    end
end

