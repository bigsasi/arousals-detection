%updateAlphaBetaOverlapping update alpha events if they overlap with beta
%events
%   [eventsAlpha] = updateAlphaBetaOverlapping(eventsAlpha, eventsBeta,
%   unitionTime) returns the updated alpha events
%
%   [eventsAlpha, updates] = updateAlphaBetaOverlapping(eventsAlpha, 
%   eventsBeta, unitionTime) returns the number of updated alpha events

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

function [eventsAlpha, updates] = updateAlphaBetaOverlapping(eventsAlpha, eventsBeta, unionTime)

updates = 0;

iniAlpha = find(diff(eventsAlpha) == 1) + 1;
endAlpha = find(diff(eventsAlpha) == -1);

% Check vector initiated by event
if endAlpha(1) < iniAlpha(1)
    endAlpha(1) = [];
end
% Check vector ended by event
if endAlpha(end) < iniAlpha(end)
    iniAlpha(end) = [];
end

for k = 1:length(iniAlpha)
   oldDuration = endAlpha(k) - iniAlpha(k);
    
   % Backwards search
    searchInterval = eventsBeta(max(iniAlpha(k)-unionTime,1):max(endAlpha(k)-1,1));

    if any(searchInterval)

        % Check the real start of the first overlapping alpha event
        iniIdx = find(searchInterval,1,'first');
        if (iniIdx == 1)
            k1 = 0;
            while eventsBeta(iniAlpha(k)-unionTime-k1-1)
                k1 = k1 + 1;
                % Check if reached beginning of recording
                if (iniAlpha(k)-unionTime-k1 == 1)
                    break;
                end
            end
            iniAlpha(k) = iniAlpha(k)-unionTime-k1;
        end
    end
        
    % Forward search
    searchInterval = eventsBeta(min(endAlpha(k)+1,length(eventsBeta)):min(endAlpha(k)+unionTime, length(eventsBeta)));

    if any(searchInterval)
        % Check the real end of the last overlapping alpha event
        endIdx = find(searchInterval,1,'last');
        if (endIdx == length(searchInterval))
            k1 = 0;
            while eventsBeta(endAlpha(k)+unionTime+k1+1)
                k1 = k1 + 1;
                % Check if reached end of recording
                if (endAlpha(k)+unionTime+k1 == length(eventsBeta))
                    break;
                end
            end
            endAlpha(k) = endAlpha(k)+unionTime+k1;
        end    
    end
    
    % Check if an update occurred
    newDuration = endAlpha(k) - iniAlpha(k);
    if ne(oldDuration, newDuration)
        % Rescore final event
        eventsAlpha(iniAlpha(k):endAlpha(k)) = 1;
        updates = updates + 1; 
    end
end