function data = detectSaccades(data,ETparams)
% Detects start and end by velocity criteria

V = data.vel;
A = data.acc;
len = length(V);

% Preallocate memory
velLabeled = bwlabel(data.velPeakIdx);
data.saccadeIdx.Idx = zeros(1,len);    % Saccade index
data.glissadeIdx.Idx = zeros(1,len);  % Glissade index

% If no saccades are detected, return
if isempty(velLabeled);
    return
end

% Process one velocity peak at the time
kk = 1;

for k = 1:max(velLabeled)
    
    %----------------------------------------------------------------------  
    % Check the saccade peak samples
    %----------------------------------------------------------------------       
    % The samples related to the current saccade
    peakIdx = find(velLabeled == k);
    
    % If the peak consists of =< minPeakSamples consequtive samples, it it probably
    % noise (1/6 or the min saccade duration)
    minPeakSamples = ceil(ETparams.minSaccadeDur/6*ETparams.samplingFreq); 
    if length(peakIdx) <= minPeakSamples, continue, end
    
    % Check whether this peak is already included in the previous saccade
    % (can be like this for glissades)
    if kk > 1
        if ~isempty(intersect(peakIdx,[find(data.saccadeIdx.Idx) find(data.glissadeIdx.Idx)]))
            continue
        end       
    end
       
    %----------------------------------------------------------------------
    % DETECT SACCADE
    %----------------------------------------------------------------------       
    
    % Detect saccade start.  AND acc <= 0
    saccadeStartIdx = find(V(peakIdx(1):-1:1) <= data.saccadeVelocityTreshold &...% vel <= global vel threshold
                                                 [diff(V(peakIdx(1):-1:1)) 0] >= 0);          % acc <= 0
    if isempty(saccadeStartIdx), continue, end
    saccadeStartIdx = peakIdx(1) - saccadeStartIdx(1) + 1;
    
    % Calculate local fixation noise (the adaptive part)
    localVelNoise = V(saccadeStartIdx:-1: max(1,ceil(saccadeStartIdx - ETparams.minFixDur*ETparams.samplingFreq)));
    localVelNoise = mean(localVelNoise) + 3*std(localVelNoise);
    localsaccadeVelocityTreshold = localVelNoise*0.3 + data.saccadeVelocityTreshold*0.7; % 30% local + 70% global
    
    % Check whether the local vel. noise exceeds the peak vel. threshold.
    if localVelNoise > data.peakDetectionThreshold, continue, end
              
    % Detect end of saccade (without glissade)
    saccadeEndIdx = find(V(peakIdx(end):end) <= localsaccadeVelocityTreshold &...             % vel <= adaptive vel threshold
                                                  [diff(V(peakIdx(end):end)) 0] >= 0);        % acc <= 0
    
    if isempty(saccadeEndIdx), continue, end      
    saccadeEndIdx = peakIdx(end) + saccadeEndIdx(1) - 1;
    
    % If the saccade contains NaN samples, continue
    if any(data.nanIdx.Idx(saccadeStartIdx:saccadeEndIdx)), continue, end
        
    % Make sure the saccade duration exceeds the minimum duration.
    saccadeLen = saccadeEndIdx - saccadeStartIdx;
    if saccadeLen/ETparams.samplingFreq < ETparams.minSaccadeDur
        continue    
    end
    
    % If all the above criteria are fulfilled, label it as a saccade.
    data.saccadeIdx.Idx(saccadeStartIdx:saccadeEndIdx) = 1;
    data.localSaccadeVelocityTreshold(kk) = localsaccadeVelocityTreshold;

        % Collect information about the saccade
    data.saccadeInfo(kk).start = saccadeStartIdx/ETparams.samplingFreq; % in ms
    data.saccadeInfo(kk).end = saccadeEndIdx/ETparams.samplingFreq; % in ms
    data.saccadeInfo(kk).duration = data.saccadeInfo(kk).end - data.saccadeInfo(kk).start;
    data.saccadeInfo(kk).amplitude = sqrt(((data.X(saccadeEndIdx)-...
                                                (data.X(saccadeStartIdx))))^2 + ...
                                               ((data.Y(saccadeEndIdx)-...
                                                (data.Y(saccadeStartIdx))))^2   );   
    data.saccadeInfo(kk).peakVelocity = max(V(saccadeStartIdx:saccadeEndIdx)); 
    data.saccadeInfo(kk).peakAcceleration = max(A(saccadeStartIdx:saccadeEndIdx)); 

    %----------------------------------------------------------------------  
    % DETECT GLISSADE (data.glissadeInfo(kk).type
    %----------------------------------------------------------------------   
    % Search only for glissade peaks in a window <= min fix duration after
    % the saccade end
    potentialGlissadeIdx = V(saccadeEndIdx: min(saccadeEndIdx + ETparams.minFixDur*ETparams.samplingFreq,len-1));
    
    % Detect glissade (low velocity criteria) 
    glissadePeakIdxW = potentialGlissadeIdx >= localsaccadeVelocityTreshold;    
    % Detect only 'complete' peaks (those with a beginning and an end)
    endIdx = find(abs(diff(glissadePeakIdxW)));   
    if length(endIdx)>1
        endIdx = endIdx(2:2:end); 
        glissadeEndWeakIdx = endIdx(end);
        nGlissadesWeak = length(endIdx);
    else
        glissadeEndWeakIdx = [];
    end
    
    % Detect glissade (high velocity criteria)
    glissadePeakIdxS = potentialGlissadeIdx >= data.peakDetectionThreshold;
    glissadeEndStrongIdx = find(glissadePeakIdxS,1,'last');       
    nGlissadesStrong = length(unique(bwlabel(glissadePeakIdxS))) - 1;
    
    % Make sure that the saccade amplitude is larger than the glissade
    % amplitued, otherwise no glissade is detected.
    if max(potentialGlissadeIdx) > max(V(saccadeStartIdx:saccadeEndIdx)) 
        glissadeEndWeakIdx = [];
        glissadeEndStrongIdx = [];
    end
    
    % If no glissade detected
   if isempty(glissadeEndWeakIdx), 
        data.glissadeInfo(kk).type = 0; 
        data.glissadeInfo(kk).duration = 0;        
   % If glissade detected    
   else
       % Detect end.
        glissadeEndIdx = saccadeEndIdx + glissadeEndWeakIdx; 
        glissadeEndIdx = glissadeEndIdx + ...
                    find(diff(V(glissadeEndIdx:end)) >= 0,1,'first') - 1; 
        glissadeIdx = saccadeEndIdx:glissadeEndIdx;
        glissadeDuration = (length(glissadeIdx))/ETparams.samplingFreq;   
    
        
        % Do not allow glissade duration > 80 ms OR
        % If the glissade contains any NaN samples, continue
        if glissadeDuration > 2*ETparams.minFixDur ||...
           any(data.nanIdx.Idx(glissadeIdx))
            data.glissadeInfo(kk).type = 0; 
            data.glissadeInfo(kk).duration = 0;
            glissadeEndStrongIdx = [];
        else
            % Collect information about the glissade
            data.glissadeInfo(kk).type = 1; % 'Weak glissade detection criteria'
            data.glissadeIdx.Idx(glissadeIdx) = 1;
            data.glissadeInfo(kk).duration = glissadeDuration;
            data.glissadeInfo(kk).numberOf = [nGlissadesWeak, nGlissadesStrong];
        end
        
   end
    
   if ~isempty(glissadeEndStrongIdx)
        data.glissadeInfo(kk).type = 2; % 'Strong glissade detection criteria'
   end                                                                                  

    kk = kk+1;
    
end

