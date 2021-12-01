function [cutoff,fixstart,fixend] = getFixations(finalweights,xpos,ypos,cutoffstd,onoffsetThresh)
% determine fixations based on finalweights from 2-means clustering

% Roy Hessels - 2014

%% input:

% finalweights              = 2-means clustering weighting
% xpos, ypos                = horizontal & vertical coordinates from ET
% cutoffstd                 = number of std above mean clustering-weight to
%                               use as fixation cutoff
% onoffsetThresh            = threshold (x*MAD of fixation) for walking 
%                               forward/back for saccade off- and onsets

%% output:

% cutoff                    = cutoff used for fixation classification
% fixstart                  = vector with fixation start indices
% fixend                    = vector with fixation end indices

%% first determine cutoff for finalweights
cutoff = mean(finalweights,'omitnan') + cutoffstd*std(finalweights,'omitnan');

% get boolean of fixations
fixbool = finalweights < cutoff;

% get indices of where fixations start and end
[fixstart,fixend] = bool2bounds(fixbool);

% for each fixation start, walk forward until recorded position is below a
% threshold of lambda*MAD away from median fixation position.
% same for each fixation end, but walk backward
for p=1:length(fixstart)
    xmedThis = median(xpos(fixstart(p):fixend(p)));
    ymedThis = median(ypos(fixstart(p):fixend(p)));
    % MAD = median(abs(x_i-median({x}))). For the 2D version, I'm using
    % median 2D distance of a point from the median fixation position. Not
    % exactly MAD, but makes more sense to me for 2D than city block,
    % especially given that we use 2D distance in our walk here
    MAD = median(hypot(xpos(fixstart(p):fixend(p))-xmedThis, ypos(fixstart(p):fixend(p))-ymedThis));
    
    thresh = MAD*onoffsetThresh;
    
    % walk until distance less than threshold away from median fixation
    % position. No walking occurs when we're already below threshold.
    i = fixstart(p);
    if i>1  % don't walk when fixation starting at start of data
        while hypot(xpos(i)-xmedThis,ypos(i)-ymedThis)>thresh
            i = i+1;
        end
        fixstart(p) = i;
    end
    
    % and now fixation end.
    i = fixend(p);
    if i<length(xpos)   % don't walk when fixation ending at end of data
        while hypot(xpos(i)-xmedThis,ypos(i)-ymedThis)>thresh
            i = i-1;
        end
        fixend(p) = i;
    end
end
