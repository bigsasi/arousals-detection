%MERGEEVENTS Merge two list of events into one
%   merged = mergeEvents(event1, event2) returns the ordered list merging
%   event1 and event2.

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

function merged = mergeEvents(event1, event2)
    nextE1 = 1;
    nextE2 = 1;
    
    length1 = length(event1);
    length2 = length(event2);
    
    if length1 == 0
        merged = event2;
        return
    end
    if length2 == 0
        merged = event1;
        return
    end
    fields1 = fieldnames(event1);
    fields2 = fieldnames(event2);
    
    merged(length1 + length2).start = 0;
    for i = 1:length1 + length2
        if nextE2 > length2 || nextE1 <= length1 && event1(nextE1).start < event2(nextE2).start
            for j = 1:length(fields1)
                merged(i).(fields1{j}) = event1(nextE1).(fields1{j});
            end
            nextE1 = nextE1 + 1;
        else
            for j = 1:length(fields2)
                merged(i).(fields2{j}) = event2(nextE2).(fields2{j});
            end
            nextE2 = nextE2 + 1;
        end
    end
end