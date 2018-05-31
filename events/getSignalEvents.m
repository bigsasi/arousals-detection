%GETSIGNALEVENTS Return events for a signal parameter given a baseline and
%certain threshold
%   [events] = getSignalEvents(signalParameter, threshold, unionTime, 
%   windowStep, baseline) returns an array with market events for a signal 
%   parameter when this paraemter is higher than a baseline multiplied by 
%   a threshold.

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

function [events] = getSignalEvents(signalParameter, threshold, unionTime, windowStep, baseline)
    numWindows = length(signalParameter);

    events = zeros(numWindows, 1, 'int8');
    events(signalParameter > (threshold * baseline)) = 1;
    % To avoid problems for detecting events starting and ending points, do not
    % allow to start or end with an event.
    events(1) = 0;
    events(end) = 0;

    % Merge events separated by less than a certain distance
    events = mergeEventsCloserThan(events, unionTime, windowStep);
end
