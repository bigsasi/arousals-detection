%% The code below is based on the methods described in the following reference(s):
% 
% [1] - I. Fernández-Varela, D. Alvarez-Estevez, E. Hernández-Pereira, V. Moret-Bonillo, 
% "A simple and robust method for the automatic scoring of EEG arousals in
% polysomnographic recordings", Computers in Biology and Medicine, vol. 87, pp. 77-86, 2017
%
% [2] - D. Alvarez-Estevez, I. Fernández-Varela, "Large-scale validation of an automatic EEG arousal detection
% algorithm using different heterogeneous databases", Sleep Medicine, vol. 57, pp. 6-14, 2019 
%
% Copyright (C) 2017-2019 Isaac Fernández-Varela
% Copyright (C) 2017-2019 Diego Alvarez-Estevez

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

function signalFilt = polymanNotchFilter(signal, sr, gain, centerFreq, bandWidth)
% Implements Notch filter from Polyman

% Check for valid settings
if ne(gain, 1)
    error('Gain should be 1 to avoid instability');
end

if ((centerFreq <= 0) || (centerFreq > sr/2))
    error('Center frequency is not in the Nyquist interval');
elseif ((bandWidth <= 0) || (bandWidth > sr/4))
    error('Not valid bandwidth');
end

r = 2 * pi * centerFreq/sr;
s = 2 * pi * (centerFreq - (bandWidth/2))/sr;
t = 1 - sqrt((cos(r)-cos(s))^2 + (sin(r)-sin(s))^2) * sqrt(exp(1)^2-1);
x = gain*cos(r);
y = gain*sin(r);

b = [gain, -2*x, x^2+y^2]; % b coeffs according to filter function description (for x(n))
a = [1, 2*t*x, -((t*x)^2 + (t*y)^2)]; % a coeffs according to filter function description (for y(n))

% We multiply a(2:end)*(-1) because filter function is implemented as a direct
% form II transposed structure
a(2:end) = -a(2:end);

signalFilt = filter(b, a, signal);
    
    