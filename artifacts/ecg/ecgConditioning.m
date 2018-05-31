%ECGCONDITIONING Removes EKG artifacts from signal
%   [SIGNAL] = ecgConditioning(SIGNAL, SRATE, EKGRATE, RELEKG, RPEAKS)
%   returns the signal without EKG artifacts. SRATE represents the signal
%   rate and EKGRATE the EKG rate. RELEKG is a vector sampled at EKGRATE
%   that indicates if the RPEAKS are reliable or not. RPEAKS is the list of
%   R-peaks inducing the artifacts. 

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

function [signalfiltered] = ecgConditioning(signal, signalRate, ecgRate, relECG, rpeaks)
%ECGCONDITIONING Conditiong of signal using ecg
%   
%   This functions aloows conditioning signal 'signal' with the ecg signal,
%   allowing to remove unwanted spikes.
%
%   signal      Signal to be conditioned
%   signalRate  Signal rate
%   ecgRate     ECG Sample Rate
%   relECG      Reliability of ECG signal
%   rpeaks      R peaks of ECG signal
    
    % Creating filtering template on the basis of the ECG peaks
    ecgArtLen = 0.5; %sec around each ECG peak
    ecgArtLenS = floor(signalRate * ecgArtLen);
    ecgArtHalf = round(ecgArtLenS / 2);
    ecgArtifact = zeros(2 * ecgArtHalf, 1);
    
    % Modify rpeaks if ecgRate ~= signalRate
    if ecgRate ~= signalRate
        rpeaksResampled = floor(rpeaks * (signalRate / ecgRate));
    else
        rpeaksResampled = rpeaks;
    end
    
    upFactor = 0.01;
    signalfiltered = signal;
    for k = 1:length(rpeaks)
        idxw = max(1, rpeaksResampled(k) - ecgArtHalf):min(length(signalfiltered), rpeaksResampled(k) + ecgArtHalf - 1);
        idxwRelECG = max(1, rpeaks(k) - ecgArtHalf):min(length(signalfiltered), rpeaks(k) + ecgArtHalf - 1);

        window = signalfiltered(idxw);
        
        if (sum(relECG(idxwRelECG)) < length(idxwRelECG))
            % It means ECG is partially not reliable on this window
            % Thus we average zeros in the ecgArtifact
            ecgArtifact = (1 - upFactor) * ecgArtifact + upFactor * zeros(size(window));

            signalfiltered(idxw) = window;
        else
            % Average current ECG artifact
            ecgArtifact = (1 - upFactor) * ecgArtifact + upFactor * window;

            % Force artifact smoothing to zero on window extremes
            taperedEcgArtifact = hann(2 * ecgArtHalf) .* ecgArtifact;

            % Remove current ECG artifact from EOG
            signalfiltered(idxw) = window - taperedEcgArtifact;
        end
    end
end