% This is a stripped and modified version of the I2MC algorithm, as used
% for Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. &
% Hessels, R.S. (2022). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.
% https://doi.org/10.3758/s13428-021-01723-1
%
% If you use this code, please cite the original paper Nyström et al.
% (2010) and the Hooge et al. (2022) paper for which the modified
% version was developed:
%
% Nyström, M., & Holmqvist, K. (2010). An adaptive algorithm for fixation,
% saccade, and glissade detection in eyetracking data. Behavior Research
% Methods, 42(1), 188–204.
%
% and
% 
% Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. & Hessels,
% R.S. (2022). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.
% https://doi.org/10.3758/s13428-021-01723-1

function episode = NH2010(dat,channel,parameters)

% get data
data.X = dat.(channel).X(:).';
data.Y = dat.(channel).Y(:).';
ETparams.samplingFreq = dat.freq;

% algo params
ETparams.blinkVelocityThreshold = 1000;             % if vel > 1000 degrees/s, it is noise or blinks
ETparams.blinkAccThreshold = 100000;                % if acc > 100000 degrees/s^2, it is noise or blinks
ETparams.peakDetectionThreshold = 100;              % Initial value of the peak detection threshold. 
ETparams.minSaccadeDur = 0; % in seconds
ETparams.minFixDur     = 0; % in seconds (used for threshold determination, 0 means any and all samples below threshold are used for determining noise of vel signal)
ETparams.sgFilterSpan  = parameters.velWindow/1000; % in seconds

% Process data
data = eventDetection(data,ETparams);

% prepare output
% In the code we find:
% data.saccadeInfo(kk).start = saccadeStartIdx/ETparams.samplingFreq;
% data.saccadeInfo(kk).end = saccadeEndIdx/ETparams.samplingFreq;
% which turns indices into time. We want to output indices, so do inverse
if ~isfield(data,'saccadeInfo')
    [episode.start, episode.end] = deal([]);
else
    episode.start = round([data.saccadeInfo.start]*ETparams.samplingFreq);
    episode.end   = round([data.saccadeInfo.end  ]*ETparams.samplingFreq);
end
