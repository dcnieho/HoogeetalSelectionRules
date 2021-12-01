function data = detectFixationNoiseLevel(data,ETparams)
 
possibleFixationIdx = ~data.InitialVelPeakIdx;
fixLabeled = bwlabel(possibleFixationIdx);

% Process one inter-peak-saccadic periods (called fixations below,
% although they are not identified as fixations yet). 
fixNoise = [];
for k = 1:max(fixLabeled)

    % The samples related to the current fixation
    fixIdx = find(fixLabeled == k);
    
    % Check that the fixation duration exceeds the minimum duration criteria. 
    if length(fixIdx)/ETparams.samplingFreq < ETparams.minFixDur
        continue    
    end
    
    % Extract the samples from the center of the fixation
    centralFixSamples = ETparams.minFixDur*ETparams.samplingFreq/6;
    fNoise = data.vel(floor(fixIdx(1)+centralFixSamples):ceil(fixIdx(end)-centralFixSamples));
    fixNoise = [fixNoise fNoise];
end

data.avgNoise = nanmean(fixNoise);
data.stdNoise = nanstd(fixNoise);

% Base the peak velocity threshold on the noise level
data.peakDetectionThreshold =  data.avgNoise + 6*data.stdNoise;
data.saccadeVelocityTreshold = data.avgNoise + 3*data.stdNoise;
data.velPeakIdx  = data.vel > data.peakDetectionThreshold;

