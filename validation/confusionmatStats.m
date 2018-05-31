%CONFUSIONMATSTATS return different performance measures
%   INPUT
%   group = true class labels
%   grouphat = predicted class labels
%
%   OR INPUT
%   stats = confusionmatStats(group);
%   group = confusion matrix from matlab function (confusionmat)
%
%   OUTPUT
%   stats is a structure array
%   stats.confusionMat
%                 Predicted Classes
%                      p'    n'
%                ___|_____|_____| 
%         Actual  p |     |     |
%        Classes  n |     |     |
%
%   Total accuracy (sum of principal diagonal divived total number of
%   cases):
%   stats.accuracy = (TP + TN)/(TP + FP + FN + TN)
%   stats.precision = TP / (TP + FP)                  % for each class label
%   stats.sensitivity = TP / (TP + FN)                % for each class label
%   stats.specificity = TN / (FP + TN)                % for each class label
%   stats.auc = (sensitivity + specificity)/2         % for each class label
%   stats.recall = sensitivity                        % for each class label
%   stats.Fscore = 2*TP /(2*TP + FP + FN)            % for each class label
%
%   TP: true positive, TN: true negative, 
%   FP: false positive, FN: false negative
% 

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

function stats = confusionmatStats(group,grouphat,classVector)

field1 = 'confusionMat';
if nargin < 2
    value1 = group;
else
    if (nargin == 3)
        value1 = confusionmat(group,grouphat,'order', classVector);
    else
        value1 = confusionmat(group,grouphat);
    end
    
end

numOfClasses = size(value1,1);
totalSamples = sum(sum(value1));
    
field2 = 'accuracy';  
value2 = sum(diag(value1))/totalSamples;

[TP,TN,FP,FN,sensitivity,specificity,precision,f_score,auc] = deal(zeros(numOfClasses,1));
for class = 1:numOfClasses
   TP(class) = value1(class,class);
   tempMat = value1;
   tempMat(:,class) = []; % remove column
   tempMat(class,:) = []; % remove row
   TN(class) = sum(sum(tempMat));
   FP(class) = sum(value1(:,class))-TP(class);
   FN(class) = sum(value1(class,:))-TP(class);
end

for class = 1:numOfClasses
    sensitivity(class) = TP(class) / (TP(class) + FN(class));
    specificity(class) = TN(class) / (FP(class) + TN(class));
    precision(class) = TP(class) / (TP(class) + FP(class));
    f_score(class) = 2*TP(class)/(2*TP(class) + FP(class) + FN(class));
    auc(class) = (sensitivity(class) + specificity(class))/2;
end

field3 = 'sensitivity';  value3 = sensitivity;
field4 = 'specificity';  value4 = specificity;
field5 = 'precision';  value5 = precision;
field6 = 'recall';  value6 = sensitivity;
field7 = 'Fscore';  value7 = f_score;
field8 = 'auc'; value8 = auc;
field9 = 'kappa'; value9 = obtainKappa(value1);
stats = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,field6,value6,field7,value7,field8,value8,field9,value9);