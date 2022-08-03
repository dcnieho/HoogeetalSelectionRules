% This is a stripped and modified version of the I2MC algorithm, as used
% for Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. &
% Hessels,  R.S. (2022). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.
% https://doi.org/10.3758/s13428-021-01723-1
%
% If you use this code, please cite the original paper Hessels et al.
% (2016) and the Hooge et al. (2022) paper for which the modified
% version was developed:
%
% Hessels, R. S., Niehorster, D. C., Kemner, C., & Hooge, I. T. C. (2016).
% Noise-robust fixation detection in eye movement data: Identification by
% two-means clustering (I2MC). Behavior Research Methods, pp 1–22.
% https://doi.org/10.3758/s13428-016-0822-1.
%
% and
% 
% Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. & Hessels,
% R.S. (2022). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.
% https://doi.org/10.3758/s13428-021-01723-1

function [fix,data,par] = I2MC(dat,channel,varargin)
% Hessels, R.S., Niehorster, D.C., Kemner, C., & Hooge, I.T.C., (2017).
% Noise-robust fixation detection in eye-movement data - Identification by 
% 2-means clustering (I2MC). Behavior Research Methods, 49(5): 1802--1823.
%
% stripped version of I2MC version v2.0.3, commit 471e49c, Feb 25, 2021

%% deal with inputs
% get inputs the user specified and throw them in the parser
if isstruct(varargin{1})
    % convert to key-value pairs
    assert(isscalar(varargin),'only one input for options is expected if options are given as a struct')
    varargin = [reshape([fieldnames(varargin{1}) struct2cell(varargin{1})].',1,[]) varargin(2:end)];
end
if isscalar(varargin) && isempty(varargin{1})
    varargin = {'freq',dat.freq};
else
    varargin = [varargin {'freq',dat.freq}];
end

data.(channel).X = dat.(channel).X;
data.(channel).Y = dat.(channel).Y;

% set defaults
% required parameters:
par.freq            = [];
% parameters with defaults:
% K-MEANS CLUSTERING
par.windowtime      = .2;       % time window (s) over which to calculate 2-means clustering (choose value so that max. 1 saccade can occur)
par.steptime        = .02;      % time window shift (s) for each iteration. Use zero for sample by sample processing
par.downsamples     = [2 5 10]; % downsample levels (can be empty)
par.downsampFilter  = 1;        % use chebychev filter when downsampling? 1: yes, 0: no. requires signal processing toolbox. is what matlab's downsampling functions do, but could cause trouble (ringing) with the hard edges in eye-movement data
par.chebyOrder      = 8;        % order of cheby1 Chebyshev downsampling filter, default is normally ok, as long as there are 25 or more samples in the window (you may have less if your data is of low sampling rate or your window is small
par.maxerrors       = 100;      % maximum number of errors allowed in k-means clustering procedure before proceeding to next file
% FIXATION DETERMINATION
par.cutoffstd       = 2;        % number of standard deviations above mean k-means weights will be used as fixation cutoff
par.onoffsetThresh  = 3;        % number of MAD away from median fixation duration. Will be used to walk forward at fixation starts and backward at fixation ends to refine their placement and stop algorithm from eating into saccades


% loop over input
checkNumeric = @(x,k) assert(isnumeric(x),'The value of ''%s'' is invalid. Expected input to be one of these types:\n\ndouble, single, uint8, uint16, uint32, uint64, int8, int16, int32, int64\n\nInstead its type was %s.',k,class(x));
checkScalar  = @(x,k) assert(isscalar(x),'The value of ''%s'' is invalid. Expected input to be a scalar.',k);
checkNumel2  = @(x,k) assert(numel(x)==2,'The value of ''%s'' is invalid. Expected input to be an array with number of elements equal to 2.',k);
checkInt     = @(x,k) assert(~any(mod(x,1)),'The value of ''%s'' is invalid. Expected input to be integer-valued.',k);
for p=1:2:length(varargin)
    key = varargin{p};
    if p+1>length(varargin)
        error('No value was given for ''%s''. Name-value pair arguments require a name followed by a value.',key);
    end
    value = varargin{p+1};
    switch key
        case {'xres','yres','freq','missingx','missingy','disttoscreen','windowtimeInterp','maxdisp','windowtime','steptime','cutoffstd','onoffsetThresh','maxMergeDist','maxMergeTime','minFixDur'}
            checkNumeric(value,key);
            checkScalar(value,key);
            par.(key) = value;
        case {'downsampFilter','chebyOrder','maxerrors','edgeSampInterp'}
            checkInt(value,key);
            checkScalar(value,key);
            par.(key) = value;
        case 'scrSz'
            checkNumeric(value,key);
            checkNumel2(value,key);
            par.(key) = value;
        case 'downsamples'
            checkInt(value,key);
            par.(key) = value;
        otherwise
            if ~ischar(key)
                error('Expected a string for the parameter name at position %d, instead the input type was ''%s''.',class(key));
            else
                error('Key "%s" not recognized',key);
            end
    end
end

% deal with required options
% if empty, user did not specify these
checkFun = @(opt,str) assert(~isempty(par.(opt)),'I2MCfunc: %s must be specified using the ''%s'' option',str,opt);
checkFun('freq', 'tracker sampling rate')

% check filter
if par.downsampFilter
    assert(exist('cheby1','file')==2,'I2MCfunc: When setting the ''downsampFilter'' option to true, the function ''cheby1'' from the signal processing toolbox is required. It appears this function is not available in your installation. Set the option to 0.')
    nSampRequired = max(1,3*par.chebyOrder)+1;  % nSampRequired = max(1,3*(nfilt-1))+1, where nfilt = chebyOrder+1
    nSampInWin    = round(par.windowtime/(1/par.freq));
    assert(nSampInWin>=nSampRequired,'I2MCfunc: Filter parameters requested with the setting ''chebyOrder'' will not work for the sampling frequency of your data. Please lower ''chebyOrder'', or set the setting ''downsampFilter'' to 0')
end
assert(~any(mod(par.freq,par.downsamples)),'I2MCfunc: Some of your downsample levels are not divisors of your sampling frequency. Change the option ''downsamples''')


%% START ALGORITHM

%% PREPARE INPUT DATA
% make sure all fields in data are columns
fs = fieldnames(data);
for f=1:length(fs)
    if isstruct(data.(fs{f}))
        fs2 = fieldnames(data.(fs{f}));
        for f2=1:length(fs2)
            data.(fs{f}).(fs2{f2}) = data.(fs{f}).(fs2{f2})(:);
        end
    else
        data.(fs{f}) = data.(fs{f})(:);
    end
end
% deal with monocular data
if isfield(data,'left')
    xpos = data.left.X;
    ypos = data.left.Y;
    missing = isnan(data.left.X) | isnan(data.left.Y);
    data.left.missing = missing;
elseif isfield(data,'right')
    xpos = data.right.X;
    ypos = data.right.Y;
    missing = isnan(data.right.X) | isnan(data.right.Y);
    data.right.missing = missing;
end

%% CALCULATE 2-MEANS CLUSTERING FOR SINGLE EYE

% get kmeans-clustering for averaged signal
[data.finalweights,stopped] = twoClusterWeighting(xpos,ypos,missing,par.downsamples,par.downsampFilter,par.chebyOrder,par.windowtime,par.steptime,par.freq,par.maxerrors);

% check whether clustering succeeded
if stopped
    fprintf('Clustering stopped after exceeding max errors\n');
    return
end
    

%% DETERMINE FIXATIONS BASED ON FINALWEIGHTS_AVG
[fix.cutoff,fix.start,fix.end] = getFixations(data.finalweights,xpos,ypos,par.cutoffstd,par.onoffsetThresh);
