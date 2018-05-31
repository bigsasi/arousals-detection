%BASELINE Mean baseline for a signal parameter
%   [baseline] = baseline(backgroundTime, windowStep, signalParameter)
%   returns a baseline for a given signal parameter stored in SIGNAL 
%   PARAMETERE using the desired background. BACKGROUNDTIME is the time in 
%   seconds for the background mean, WINDOWSTEP is the displacement between 
%   two consecutive measures in SIGNALPARAMETER. 

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

function [baseline] = baseline(backgroundTime, windowStep, signalParameter)
    backgroundSize = round(backgroundTime / windowStep);
	backgroundWindows = length(1:windowStep:backgroundTime);
    numWindows = length(signalParameter);
    
    rawBaseline = zeros(backgroundSize, numWindows);
    background(1:backgroundSize) = mean(signalParameter(1:backgroundWindows));
    next = 1;
    for k = 1:numWindows
        rawBaseline(:, k) = background;
        next = next + 1;
        background(next) = signalParameter(k);
        if next == backgroundSize, next = 0; end
    end
    baseline = sum(rawBaseline, 1)' / backgroundSize;
end

