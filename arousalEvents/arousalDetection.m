%arousalDetection Find arousals using the eeg and emg signal and the
%hypnogram
%   [AR] = arousalDetection(EEG, EMG) returns the list of
%   arousals found in the input eeg, using the emg to improve the detection.
%
%   [AR] = arousalDetection(EEG, EMG, HYPNOGRAM) also uses the hypnogram to
%   improve the detection.
%
%   [AR, F] = arousalDetection(EEG, EMG, HYPNOGRAM) also returns the list
%   of potential arousals with the info of why they are discarded.
%
%   [AR, F, S] = arousalDetection(EEG, EMG, HYPNOGRAM) also returns the
%   list of spindles found in the discarded arousals. 
%
%   [AR, F, S, EEG] = arousalDetection(EEG, EMG, HYPNOGRAM) also return the
%   EEG signal including the power events used to detect the arousals. 
%
%   [AR, F, S, EEG, EMG] = arousalDetection(EEG, EMG, HYPNOGRAM) also
%   return the EMG signal including the amplitude events used to improve
%   the detection. 

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

function [arousals, removedArousals, spindles, eeg, emg] = arousalDetection(eeg, emg, hypnogram)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    logger = logging.Logger.getLogger('ArousalDetection');
    
    if nargin < 3
        logger.info('Performance can improve using an expert hypnogram');
        hypnogram = ones(ceil(length(eeg.raw) / eeg.rate / 30), 1);
    end

    logger.info('Finding power events in EEG');
    eeg = findPowerEvents(eeg);

    logger.info('Finding amplitude events in EMG');
    emg = findAmplitudeEvents(emg);
    emg.amplitudeUpdate = signalAmplitude(emg.raw, emg.rate, Configuration.instance().updateWithEmgWindowSize, Configuration.instance().updateWithEmgWindowStep);
    
    [detectedArousals, removedArousals, spindles] = findArousalEvents(eeg, emg, hypnogram);
    
    arousals = detectedArousals;
end