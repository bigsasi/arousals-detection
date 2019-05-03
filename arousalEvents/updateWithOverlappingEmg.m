%updateWithOverlappingEmg Update the event using the emg
%signal. 
%   EVENT = updateWithOverlappingEmg(EVENT, EEG) returns the 
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

function event = updateWithOverlappingEmg(event, emg)
%     logger = logging.Logger.getLogger('updateWithOverlappingEmg');

    windowSize = Configuration.instance().updateWithEmgWindowSize;
    windowStep = Configuration.instance().updateWithEmgWindowStep;
    ampRate = 1 / windowStep;
    threshold = Configuration.instance().updateWithEmgThreshold;
    
    deltaSamples = 15 * emg.rate;
    leftThresholdSample = event.startSample - deltaSamples;
    rigthThresholdSample = event.endSample + deltaSamples;

    aux = signalAmplitude(emg.raw(leftThresholdSample:rigthThresholdSample), emg.rate, windowSize, windowStep);
    aux = sort(aux);

    % Select a threshold using the P50 amplitude
    amplitudeThreshold = threshold * aux(round(0.5 * length(aux)));        

    startSample = round(event.start * ampRate); 
    endSample = round((event.start + event.duration) * ampRate);

    if sum(emg.amplitudeUpdate(endSample - 10:endSample) > amplitudeThreshold) >= 8
        for emgEventEndSample = endSample:length(emg.amplitudeUpdate)
            % Finish the update if amplitude is lower than the threshold
            if emg.amplitudeUpdate(emgEventEndSample) <= amplitudeThreshold
                break;
            end
        end
    else
        emgEventEndSample = endSample;
    end
    
    emgEventLength = emgEventEndSample - startSample;
    duration = emgEventLength / ampRate;

    if duration > event.duration + 0.1
        event.info{end + 1} = ['EMG update +' num2str(duration - event.duration) ' s'];
        event.endSample = event.startSample + round(duration * emg.rate);
        event.duration = duration;    
    end    
end