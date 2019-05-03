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

function out = ExpInteger(input, Y0, A)

if (input > 0)
    out = Y0 * exp(A * input);
elseif (input < 0)
    out = (-1)*Y0 * exp((-1) * A * input);
else
    out = 0;
end