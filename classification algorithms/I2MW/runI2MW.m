clear all; close all; clc

% set parameters
params.windowDuration       = 21;       % ms
params.amplitudeThreshold   = 0.1;      % deg
params.allowMissing         = false;    % if true, allows missing during saccades

% import data
outputMatrix    = readNumericFile('testdata.txt',5,1);
data.t          = outputMatrix(:,1);    % time signal
data.left.X     = outputMatrix(:,2);    % horizontal gaze signal of left eye
data.left.Y     = outputMatrix(:,3);    % vertical gaze signal of left eye
data.right.X    = outputMatrix(:,4);    % horizontal gaze signal of right eye
data.right.Y    = outputMatrix(:,5);    % vertical gaze signal of right eye
% add other info about data
data.freq       = 1000;                 % sampling frequency (Hz)

% do saccade classification
% the output is a struct containing the starts and ends of saccades
% (sample numbers)
episodes        = I2MW(data,'left',params);

% convert from sample number to time (ms)
episodes.startT = data.t(episodes.start);
episodes.endT   = data.t(episodes.end);

% print saccade start and end on command window
fprintf('%d saccades classified\n',length(episodes.startT));
fprintf('saccade start -- end: duration (all in ms)\n');
for p=1:length(episodes.startT)
    fprintf('%3d: %4.0f -- %4.0f: %4.0f\n',p,episodes.startT(p),episodes.endT(p),episodes.endT(p)-episodes.startT(p));
end