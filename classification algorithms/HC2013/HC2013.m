% This is a stripped and modified version of the HC2013 algorithm, as used
% for Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. &
% Hessels, R.S. (in press). Fixation classification: how to merge and
% select fixation candidates. Behavior Research Methods.
%
% If you use this code, please cite the original paper Hooge & Camps (2013)
% and the Hooge et al. (in press) paper for which the modified version was
% developed:
%
% Hooge, I. T. C., & Camps, G. (2013). Scan path entropy and arrow plots:
% Capturing scanning behavior of multiple observers. Frontiers in
% Psychology, 4, 996. https://doi.org/10.3389/fpsyg.2013.00996
%
% and
% 
% Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. & Hessels,
% R.S. (in press). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.

function [fix] = HC2013(data,channel,parameters)

time = data.t;
xd = data.(channel).X;
yd = data.(channel).Y;


% algorithm parameters
f       = parameters;
f.freq  = data.freq;
f.tc    = 1000/f.freq;          % sample interval
span    = parameters.velWindow/1000*f.freq;
assert(mod(span,1)==0,'window duration %.1f ms is not multiple of sample interval %.1f ms',parameters.velWindow,1000/f.freq)
assert(mod(span,2)==1,'Window duration %.1f ms should correspond to an uneven number of samples, not %d samples',parameters.velWindow, span);
f.dn    = (span-1)/2;
fprintf('Window duration used: %1$.1f ms, %2$d samples (%3$d+1+%3$d)\n',parameters.velWindow,span,f.dn);

% get velocity
vel1    = getvelacc(xd,f.dn,f.freq);
vel2    = getvelacc(yd,f.dn,f.freq);
mvel    = sqrt(vel1.^2 + vel2.^2);

% prep output
fix.start   = [];
fix.end     = [];

% run algo
if ~isempty(mvel)
    mvel                        = [mvel(f.dn+1)*ones(f.dn,1);mvel(1+f.dn:end-f.dn);mvel(end-f.dn)*ones(f.dn,1)];

    [fmark,thr2,meanvel,stdvel]	= detectfixaties2020_DN(mvel,f);
    
    % make output
    fix.start  = fmark(1:2:end);
    fix.end    = fmark(2:2:end);
end

