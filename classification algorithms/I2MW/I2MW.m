% This is the I2MW algorithm, as developed for and used in Hooge, I.T.C.,
% Niehorster, D.C., Nyström, M., Andersson, R. & Hessels, R.S. (in press).
% Fixation classification: how to merge and select fixation candidates.
% Behavior Research Methods.
%
% If you use this code, please cite the Hooge et al. (in press) paper for
% which this algorithm was developed:
% 
% Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. & Hessels,
% R.S. (in press). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.

function [sac] = I2MW(data,channel,parameters)
% this algorithm deploys two moving windows at a constant separation over
% the data. If the difference in mean position between the two windows
% exceeds a threshold, all samples between the two windows are marked as
% part of a saccade.

xd = data.(channel).X(:);
yd = data.(channel).Y(:);


% algorithm parameters, convert from ms to samples, where needed
span    = parameters.windowDuration/1000*data.freq;
assert(mod(span,1)==0,'window duration %.1f ms is not multiple of sample interval %.1f ms',parameters.windowDuration,1000/data.freq)
assert(mod(span,2)==1,'Window duration %.1f ms should correspond to an uneven number of samples, not %d samples',parameters.windowDuration, span);
windowLength        = (span-1)/2;
windowSeparatation  = 1;
fprintf('Window duration used: %1$.1f ms, %2$d samples (%3$d+%4$d+%3$d)\n',parameters.windowDuration,span,windowLength,windowSeparatation);
amplitudeThreshold  = parameters.amplitudeThreshold;
allowMissing        = parameters.allowMissing;


% figure out window positions
nSamp = length(xd);
lastIdx = nSamp - 2*windowLength-windowSeparatation + 1;

w1 = [1:windowLength];
w2 = w1+windowLength+windowSeparatation;
sep= windowLength+[1:windowSeparatation];
isSac = false(size(xd));
for p=1:lastIdx
    xw1 = xd(p-1+w1);
    yw1 = yd(p-1+w1);
    xw2 = xd(p-1+w2);
    yw2 = yd(p-1+w2);
    
    p1 = median([xw1 yw1],1,'omitnan');
    p2 = median([xw2 yw2],1,'omitnan');
    
    if hypot(p1(1)-p2(1),p1(2)-p2(2)) > amplitudeThreshold
        % mark samples in between the windows as saccade
        isSac(p-1+sep) = true;
    end
end

% turn boolean into starts and ends
[sac.start, sac.end] = bool2bounds(isSac);
% remove episodes containing missing, if wanted
if ~allowMissing
    for p=length(sac.start):-1:1
        if any(isnan(xd(sac.start(p):sac.end(p))))
            sac.start(p)= [];
            sac.end(p)  = [];
        end
    end
end