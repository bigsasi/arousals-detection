%CONFIGURATION class for loading the configuration parameters.

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

classdef Configuration
    %CONFIGURATION Summary of this class goes here
    %   Detailed explanation goes here
    properties (Access = private)
        singleton
    end
    properties
        powerEventsAlphaThreshold
        powerEventsAlphaUnion
        powerEventsBetaThreshold 
        powerEventsBetaUnion
        powerEventsWindowStep
        powerEventsWindowSize
        powerEventsBackgroundTime
        
        amplitudeEventsThreshold
        amplitudeEventsWindowStep
        amplitudeEventsWindowSize
        amplitudeEventsBackgroundTime
        amplitudeEventsMergeNear
        
        updateWithAlphaEvents
        updateWithAlphaEventsUnionTime 
        updateWithBetaEvents 
        updateWithBetaEventsUnionTime 
        
        updateWithPowerBetaThreshold 
        updateWithPowerBetaThreshold2 
        updateWithPowerBetaLimitTime 
        updateWithPowerBetaIncreaseTime 
        updateWithPowerBetaPreviousSeconds 
        
        updateWithAmplitudeThreshold 
        
        updateWithFutureEmgLimit 
        updateWithFutureEmgThreshold 
        
        updateWithPowerAlphaThreshold 
        updateWithPowerAlphaThreshold2 
        updateWithPowerAlphaLimitTime 
        updateWithPowerAlphaIncreaseTime 
        updateWithPowerAlphaPreviousSeconds 
        
        updateWithEmgWindowSize 
        updateWithEmgWindowStep 
        updateWithEmgThreshold 
    end
    
    methods(Access=private)
        function newObj = Configuration()
            ini = IniConfig();
            ini.ReadFile('config.ini');
            newObj.powerEventsAlphaThreshold = ini.GetValues('PowerEvents', 'AlphaThreshold');
            newObj.powerEventsAlphaUnion = ini.GetValues('PowerEvents', 'AlphaUnion');
            newObj.powerEventsBetaThreshold = ini.GetValues('PowerEvents', 'BetaThreshold');
            newObj.powerEventsBetaUnion = ini.GetValues('PowerEvents', 'BetaUnion');
            newObj.powerEventsWindowStep = ini.GetValues('PowerEvents', 'WindowStep');
            newObj.powerEventsWindowSize = ini.GetValues('PowerEvents', 'WindowSize');
            newObj.powerEventsBackgroundTime = ini.GetValues('PowerEvents', 'BackgroundTime');

            newObj.amplitudeEventsThreshold = ini.GetValues('AmplitudeEvents', 'Threshold');
            newObj.amplitudeEventsWindowStep = ini.GetValues('AmplitudeEvents', 'WindowStep');
            newObj.amplitudeEventsWindowSize = ini.GetValues('AmplitudeEvents', 'WindowSize');
            newObj.amplitudeEventsBackgroundTime = ini.GetValues('AmplitudeEvents', 'BackgroundTime');
            newObj.amplitudeEventsMergeNear = ini.GetValues('AmplitudeEvents', 'MergeNear');

            newObj.updateWithAlphaEvents = ini.GetValues('UpdateWith', 'AlphaEvents');
            newObj.updateWithAlphaEventsUnionTime = ini.GetValues('UpdateWith', 'AlphaEventsUnionTime');
            newObj.updateWithBetaEvents = ini.GetValues('UpdateWith', 'BetaEvents');
            newObj.updateWithBetaEventsUnionTime = ini.GetValues('UpdateWith', 'BetaEventsUnionTime');

            newObj.updateWithPowerBetaThreshold = ini.GetValues('UpdateWithPower', 'BetaThreshold');
            newObj.updateWithPowerBetaLimitTime = ini.GetValues('UpdateWithPower', 'BetaLimitTime');
            newObj.updateWithPowerBetaIncreaseTime = ini.GetValues('UpdateWithPower', 'BetaIncreaseTime');
            newObj.updateWithPowerBetaPreviousSeconds = ini.GetValues('UpdateWithPower', 'BetaPreviousSeconds');

            newObj.updateWithAmplitudeThreshold = ini.GetValues('UpdateWithAmplitude', 'Threshold');

            newObj.updateWithFutureEmgLimit = ini.GetValues('UpdateWithFutureEmg', 'Limit');
            newObj.updateWithFutureEmgThreshold = ini.GetValues('UpdateWithFutureEmg', 'Threshold');

            newObj.updateWithPowerAlphaThreshold = ini.GetValues('UpdateWithPower', 'AlphaThreshold');
            newObj.updateWithPowerAlphaLimitTime = ini.GetValues('UpdateWithPower', 'AlphaLimitTime');
            newObj.updateWithPowerAlphaIncreaseTime = ini.GetValues('UpdateWithPower', 'AlphaIncreaseTime');
            newObj.updateWithPowerAlphaPreviousSeconds = ini.GetValues('UpdateWithPower', 'AlphaPreviousSeconds');

            newObj.updateWithEmgWindowSize = ini.GetValues('UpdateWithEmg', 'WindowSize');
            newObj.updateWithEmgWindowStep = ini.GetValues('UpdateWithEmg', 'WindowStep');
            newObj.updateWithEmgThreshold = ini.GetValues('UpdateWithEmg', 'Threshold');
        end
    end
   
    methods(Static)
        function obj = instance()
            persistent singleton
            if isempty(singleton)
                obj = Configuration();
                singleton = obj;
            else
                obj = singleton;
            end
        end
   end
    
end

