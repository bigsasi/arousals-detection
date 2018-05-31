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

function signalFilt = polymanHighPassFilter(signal, sr, gain, hpFreq)
% Implements high pass filter from Polyman

% Check for valid settings
if ((hpFreq <= 0) || (hpFreq > sr/2))
    disp('Error');
    return;
end

r = 2 * pi * hpFreq * ((1/sr)/2);
t = gain / (r + 1);

b = [t, - t]; % b coeffs according to filter function description (for x(n))
a = [1, -(r-1)/(r+1)]; % a coeffs according to filter function description (for y(n))

% We multiply a(2:end)*(-1) because filter function is implemented as a direct
% form II transposed structure
a(2:end) = -a(2:end);

signalFilt = filter(b, a, signal);
    