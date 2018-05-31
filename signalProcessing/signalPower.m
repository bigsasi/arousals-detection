%SIGNALPOWER study the power of a signal
%   [POWERS] = signalPower(SIGNAL, WINDOWSIZE, SIGNALRATE, ...
%   WINDOWSTEP, LCUTOFF, HCUTOFF) returns a list of vectors containing the
%   power values for SIGNAL, obtained with a sliding window of WINDOWSIZE
%   seconds with step WINDOWSTEP. Each power band is specified with the
%   lower limit in LCUTOFF and the upper limit y HCUTOFF.

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

function [varargout] = signalPower(signal, windowSize, signalRate, ...
    windowStep, lcutoff, hcutoff)
%SIGNALPOWER This function obtains the power for each window of the input
%signal, given a window size and window step. Can obtain the power in multiple
%bands at once

    blockSize = 10000; % In samples, controls size of block processing for FFT (memory/speed compromise)

    windowStepSamples = signalRate * windowStep;
    windowSizeSamples = windowSize * signalRate;
    signalLength = length(signal);
    longitudLowFrecs = length(lcutoff);
    hammingWindow = hamming(windowSizeSamples);

    if (any(hcutoff > (signalRate / 2)) || any(lcutoff > (signalRate / 2)))
        fprintf(1, '\nError: Sample rate is %.2f Hz. => only aplicable in [0 - %.2f ]Hz', signalRate, signalRate/2);
        return;
    end

    if ne(longitudLowFrecs, length(hcutoff))
        fprintf(1, '\nError when applying the band-pass filter: The number of the lower cutting points differs from the number of upper cutting points');
        return;
    end
    
    numWindows = length(1:windowStepSamples:signalLength);
    
    % Calculate the number of unique points (depending on nfft odd or even)
    NumUniquePts = ceil((windowSizeSamples+1)/2);
    % We construct the evenly spaced frequency vector with NumUniquePts to
    % represent x-axis
    step = signalRate / windowSizeSamples;
    f = 0:step:(NumUniquePts - 1) * step;
    % Index values that are equal or higher to corresponding frequency limits
    f1 = ceil(lcutoff / step) + 1;
    f2 = ceil(hcutoff / step) + 1;
    f2(f2 > length(f)) = length(f);

    varargout{nargout} = 0; 
    
    windowsLimit = length(1:windowStepSamples:signalLength - windowSizeSamples);
    currentSample = 1;
    for i = 1:blockSize:windowsLimit
        windowsSize = min(blockSize, windowsLimit - i + 1);
        windows = zeros(windowSizeSamples, windowsSize);    
        for currentWindow = 1:windowsSize
            sample = floor(currentSample);
            limit = sample + windowSizeSamples - 1;
            windows(:, currentWindow) = signal(sample:limit) .* hammingWindow;
            currentSample = currentSample + windowStepSamples;
        end

        fftx(:, i:i + windowsSize - 1) = fft(windows, windowSizeSamples);
    end
    fftx(:, windowsLimit + 1:numWindows) = 0;
    currentWindow = windowsLimit;
    
    % Since FFT is symmetric, throw away second half
    mx = abs(fftx(1:NumUniquePts, :));
    % Take the square magnitude, thus to obtain the power (without scaling)
    mx = mx .^ 2;
    % Since we dropped half the FFT, we multiply mx by 2 to keep the same energy. 
    % The DC component (fftx(1)) and Nyquist component (fftx(1+nfft/2), if it
    % exists) are unique and should not be multiplied by 2.
    limit = NumUniquePts - 1;
    % odd nfft excludes Nyquist point
    if rem(windowSizeSamples, 2), limit = NumUniquePts; end
    
    mx = mx / windowSizeSamples;
    
    for k = 1:nargout
        if f1(k) > 1 && f2(k) <= limit
            varargout{k} = sum(mx(f1(k):f2(k), :), 1)' * 2;     
        elseif f1(k) > 1 && f2(k) > limit
            varargout{k} = sum(mx(f1(k):limit, :), 1)' * 2 + sum(mx(limit + 1:f2(k), :), 1)';
        elseif f2(k) < limit
            varargout{k} = sum(mx(1, :), 1)' + sum(mx(2:f2(k), :), 1)' * 2;
        else
            varargout{k} = sum(mx(1, :), 1)' + sum(mx(2:limit, :), 1)' * 2 + sum(mx(limit + 1:f2(k), :), 1)';
        end
        varargout{k} = varargout{k} / windowSizeSamples;
        varargout{k}(currentWindow + 1:numWindows) = varargout{k}(currentWindow);
    end
end