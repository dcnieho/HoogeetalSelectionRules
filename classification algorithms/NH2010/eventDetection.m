function data = eventDetection(data,ETparams)

% Calculate velocity and acceleration
%-------------------------------------
data = calVelAcc_sgolay(data,ETparams);

% Detect blinks and noise
%-------------------------------------
data = detectAndRemoveNoise(data,ETparams);

% iteratively find the optimal noise threshold
%-------------------------------------
data.peakDetectionThreshold = ETparams.peakDetectionThreshold;
oldPeakT = inf;
while abs(data.peakDetectionThreshold -  oldPeakT) > 1
    
    oldPeakT  = data.peakDetectionThreshold;
    
    % Detect peaks in velocity (> X degrees/second)
    data = detectVelocityPeaks(data);
    
    % Find fixation noise level (0.7*global fixation noise +
    % 0.3*local fixation)
    data = detectFixationNoiseLevel(data,ETparams);
    
end

% Detect saccades (with peak detection threshold (v < v_avg_noise + 3*v_std_noise))
% and glissades
%-------------------------------------
data = detectSaccades(data,ETparams);