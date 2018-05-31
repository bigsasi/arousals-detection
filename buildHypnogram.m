%BUILDHYPNOGRAM returns an hypnogram from the sleep stages annotations
%   [hypnogram] = buildHypnogram(annotations, duration) returns a vector
%   representing the hypnogram from the list of annotations. 

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

function [hypnogram] = buildHypnogram(annotations, duration)
%BUILDHYPNOGRAM Summary of this function goes here
%   Detailed explanation goes here
    logger = logging.Logger.getLogger('BuildHypnogram');
    numEpochs = ceil(duration / 30);

    hypnogram = zeros(1, numEpochs);
    hypnogram(:) = -1;
    
    for i = 1:length(annotations.SleepstageW) 
        event = annotations.SleepstageW(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 0;
    end
    for i = 1:length(annotations.SleepstageN1)
        event = annotations.SleepstageN1(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 1;
    end
    for i = 1:length(annotations.Sleepstage1)
        event = annotations.Sleepstage1(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 1;
    end
    for i = 1:length(annotations.SleepstageN2) 
        event = annotations.SleepstageN2(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 2;
    end
    for i = 1:length(annotations.Sleepstage2) 
        event = annotations.Sleepstage2(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 2;
    end
    for i = 1:length(annotations.SleepstageN3) 
        event = annotations.SleepstageN3(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 3;
    end
    for i = 1:length(annotations.Sleepstage3) 
        event = annotations.Sleepstage3(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 3;
    end
    for i = 1:length(annotations.Sleepstage4) 
        event = annotations.Sleepstage4(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 3;
    end
    for i = 1:length(annotations.SleepstageR) 
        event = annotations.SleepstageR(i);
        epoch = floor(event.start / 30) + 1;
        total = floor(event.duration / 30);
        hypnogram(epoch:epoch + total - 1) = 5;
    end
    
    if any(hypnogram == -1) 
        logger.warning('No sleep stage for some epochs');
    end
end

