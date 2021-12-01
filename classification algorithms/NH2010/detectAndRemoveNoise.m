function data = detectAndRemoveNoise(data,ETparams)
% Detects and removes un-physiological movement (which derives from noise
% and blinks)

data.nanIdx.Idx = zeros(1,length(data.X));

V = data.vel;
V_threshold = median(data.vel)*2;

% Detect possible blinks and noise (where samples are nan or if the eyes
% move too fast)
blinkIdx = isnan(data.X) |...
        data.vel > ETparams.blinkVelocityThreshold |...
        abs(data.acc) > ETparams.blinkAccThreshold;

% Set possible blink and noise index to '1'
data.nanIdx.Idx(blinkIdx) = 1;

% Label blinks or noise
blinkLabeled = bwlabel(blinkIdx);

% Process one blink or noise period at the time
for k = 1:max(blinkLabeled)

    % The samples related to the current event
    b = find(blinkLabeled == k);
      
    % Go back in time to see where the blink (noise) started
    sEventIdx = find(V(b(1):-1:1) <= V_threshold);
    if isempty(sEventIdx), continue, end
    sEventIdx = b(1) - sEventIdx(1) + 1;
    data.nanIdx.Idx(sEventIdx:b(1)) = 1;      
    
    % Go forward in time to see where the blink (noise) started    
    eEventIdx = find(V(b(end):end) <= V_threshold);
    if isempty(eEventIdx), continue, end    
    eEventIdx = (b(end) + eEventIdx(1) - 1);
    data.nanIdx.Idx(b(end):eEventIdx) = 1;
    
end

temp_idx = find(data.nanIdx.Idx);
if temp_idx/length(V) > 0.20
    disp('Warning: This trial contains > 20 % noise+blinks samples')
    data.NoiseTrial = 0;
else
    data.NoiseTrial = 1;
end
data.vel(temp_idx) = nan;
data.acc(temp_idx) = nan;
data.X(temp_idx) = nan;
data.Y(temp_idx) = nan;


