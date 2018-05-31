%ECGREALIABILITYCHECK Returns a vector indicating the EKG signal
%reliability
%   REL = ecgReliabilityCheck(ECG, SR) returns a vector sampled at SR
%   indicating if the ECG (also sampled at SR) is reliable or not.

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

function relECG = ecgReliabilityCheck(ecg, sr)

useAmplitudeNormalization = 1;

relECG = ones(size(ecg));

%% Notch filter

fsignal = polymanNotchFilter(ecg, sr, 1, 50, 1);

%% Band-pass filter
fsignal = polymanBandPassFilter(fsignal, sr, 1, 8, 16);

%% Derivative filter
fsignal = diff(fsignal);

%% Rectifier (abs)
fsignal = abs(fsignal);

%% Mean filter
n = round(0.08 * sr);
meanFilt = 1/n * ones(1, n);
fsignal = conv(fsignal, meanFilt);
% Avoid abrupt starts/endings due to zero padding
fsignal(1:n) = fsignal(n);
fsignal(end-n+1:end) = fsignal(end-n+1);
% Compensation for filter delay
fsignal = fsignal(floor(n/2):(length(fsignal)-floor(n/2)));

%% Amplitude normalization
n = round(0.5 * sr);
meanFilt = 1/n * ones(1, n);
fnorm = conv(fsignal, meanFilt);
% Avoid abrupt starts/endings due to zero padding
fnorm(1:n) = fnorm(n);
fnorm(end-n+1:end) = fnorm(end-n+1);
% Compensation for filter delay
fnorm = fnorm(floor(n/2):(length(fnorm)-floor(n/2)));

% Workaround to fix possible duration discrepancies due to rounding
if length(fnorm) < length(ecg)
    fnorm = [fnorm; ones(length(ecg)-length(fnorm), 1)*fnorm(end)];
end
if length(fsignal) < length(ecg)
    fsignal = [fsignal; ones(length(ecg)-length(fsignal), 1)*fsignal(end)];
end

if useAmplitudeNormalization
    % Actual normalization
    fsignal = fsignal ./ fnorm;
end

%% Peak detection

wsize = 10; %seconds
step = 10; % In seconds
steps = step * sr; % algorithm step in samples
wsizes = wsize * sr; % window size in samples
numWindows = floor((length(fsignal) - wsizes)/steps);

peakRegions = zeros(size(fsignal));
valleyRegions = zeros(size(fsignal));
baseline = zeros(size(fsignal));
for k = 1:numWindows
    idxw = (k -1)*steps + 1:(k-1)*steps + wsizes;
    window = fsignal(idxw);
        
    baseline(idxw) = mean(window);
    peakRegions(idxw) = (window > baseline(idxw));
    valleyRegions(idxw) = (window <= baseline(idxw));
end
inisP = find(diff(peakRegions) == 1) + 1;
endsP = find(diff(peakRegions) == -1);

if peakRegions(1)
    inisP = [1; inisP];
end
if peakRegions(end)
    endsP = [endsP; length(peakRegions)];
end

if length(inisP) > length(endsP)
    inisP = inisP(1:end-1);
end
if length(endsP) > length(inisP)
    endsP = endsP(2:end);
end

peaks = zeros(length(inisP),1);
idxPeaks = zeros(length(inisP),1);
signalPeaks = nan(size(fsignal));
for k = 1:length(inisP)
    [peaks(k), idxPeak] = max(fsignal(inisP(k):endsP(k)));
    idxPeaks(k) = inisP(k) + idxPeak;
    signalPeaks(idxPeaks(k)) = fsignal(idxPeaks(k));
end

valleys = zeros(length(inisP) - 1,1);
idxValleys = zeros(length(inisP) - 1,1);
signalValleys = nan(size(fsignal));
for k = 2:length(inisP)
    [valleys(k-1), idxValley] = min(fsignal(endsP(k-1):inisP(k)));
    idxValleys(k-1) = endsP(k-1) + idxValley;
    signalValleys(idxValleys(k-1)) = fsignal(idxValleys(k-1));
end

%% Histogram analysis
wsize = 30; %seconds
step = 30; % In seconds
steps = step * sr; % algorithm step in samples
wsizes = wsize * sr; % window size in samples
numWindows = floor((length(fsignal) - wsizes)/steps);

data.idxPeaks = idxPeaks;
data.peaks = peaks;

for k = 1:numWindows
    startSec = (k -1)*step + 1;
    idxw = (k -1)*steps + 1:(k-1)*steps + wsizes;
    
    relECG(idxw) = analyzeReliability(startSec, wsize, sr, data);
end

function [reliable] = analyzeReliability(testsec, numSecs, sr, data)

testoffset = testsec * sr;
offlength = sr * numSecs;

mPeaks = data.peaks(data.idxPeaks >= testoffset & data.idxPeaks <= testoffset+offlength-1);

% Attempting classification

% We assume we perform an analysis on a analysisBPM second basis, then we assume we
% can trust the analysisBPM higher peaks. We are thus assuming an equivalent to analysisBPM bpm.
% This is not risky, if we have higher rate, then still the analysisBPM first peaks
% should be ok (should have little variability). If we have less than
% analysisBPM
% peaks in a minute basis (bpm < analysisBPM), then we still can reason with the peaks we have.
% The reason to reason on a "subset" of the total number of peaks is that
% sometimes P waves are included (it depends on the accuracy of the
% threshold), which generates two clusters of peaks (R waves, P waves) and
% this can hamper analysis of variability (P waves are introducing
% artificial variability in the R waves)
analysisBPM = 40; %bpm
peaks2pick = round(analysisBPM * numSecs/60);
temp = sort(mPeaks, 1, 'descend');
targetPeaks = temp(1:min(length(mPeaks), peaks2pick));
analysisRatio = std(targetPeaks)/mean(targetPeaks);

reliabilityThreshold = 0.05; % Empirically determined
if (analysisRatio > reliabilityThreshold)
    reliable = 0;
else
    reliable = 1;
end
