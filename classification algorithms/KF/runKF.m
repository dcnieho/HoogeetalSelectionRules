clear all; close all; clc

% set parameters
params.chi_threshold        = 5;
params.deviation            = 3000;
params.window_size          = 2;    % # samples

% import data
outputMatrix    = readNumericFile('testdata.txt',5,1);
data.t          = outputMatrix(:,1);    % time signal
data.left.X     = outputMatrix(:,2);    % horizontal gaze signal of left eye
data.left.Y     = outputMatrix(:,3);    % vertical gaze signal of left eye
data.right.X    = outputMatrix(:,4);    % horizontal gaze signal of right eye
data.right.Y    = outputMatrix(:,5);    % vertical gaze signal of right eye
% add other info about data
data.freq       = 1000;                 % sampling frequency (Hz)

% do fixation classification
% the output is a struct containing the starts and ends of fixations
% (sample numbers)
episodes        = KF(data,'left',params,1);

% convert from sample number to time (ms)
episodes.startT = data.t(episodes.start);
episodes.endT   = data.t(episodes.end);

% print fixation start and end on command window
fprintf('%d fixations classified\n',length(episodes.startT));
fprintf('fixation start -- end: duration (all in ms)\n');
for p=1:length(episodes.startT)
    fprintf('%3d: %4.0f -- %4.0f: %4.0f\n',p,episodes.startT(p),episodes.endT(p),episodes.endT(p)-episodes.startT(p));
end