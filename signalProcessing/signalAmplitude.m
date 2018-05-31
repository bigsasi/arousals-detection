%SIGNALAMPLITUDE study the amplitude of a signal
%   [AMPLITUDES] = signalAmplitude(SIGNAL, SIGNALRATE, WINDOWSIZE, 
%   WINDOWSTEP) returns a vector with the amplitudes of SIGNAL (sampled at
%   SIGNALRATE) using sliding windows of WINDOWSIZE seconds with a shifting
%   step of WINDOWSTEP seconds. 

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

function [amplitudes] = signalAmplitude(signal, signalRate, windowSize, windowStep)
    windowStepSamples = signalRate * windowStep;
    windowSizeSamples = signalRate * windowSize;
    signalLength = length(signal);
    
    completeWindows = floor((signalLength - windowSizeSamples) / windowStepSamples);
    
    % Indexes to select signal segments based on windowSize and windowStep
    windowsIdx = repmat(1:floor(windowSizeSamples), [completeWindows 1]);
    windowsIdx = windowsIdx + repmat(floor((0:completeWindows - 1)' * windowStepSamples), [1 floor(windowSizeSamples)]);

    % Select the segments and compute amplitude
    selectedWindows = signal(windowsIdx);
    amplitudes = max(signal(windowsIdx), [], 2) - min(selectedWindows, [], 2);
    
    currentWindow = completeWindows;
    currentSample = currentWindow * windowStepSamples + 1;
    
    for currentSample = currentSample:windowStepSamples:signalLength
        currentWindow = currentWindow + 1;
        sample = round(currentSample);
        limit = signalLength;
        
        selectedWindow = signal(sample:limit);
        
        amplitudes(currentWindow) = max(selectedWindow) - min(selectedWindow);
    end
end

