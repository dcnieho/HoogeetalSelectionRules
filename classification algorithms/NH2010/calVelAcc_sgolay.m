function data = calVelAcc_sgolay(data,ETparams)
% Lowpass filter window length
smoothInt = ETparams.sgFilterSpan; % in seconds

% Span of filter
span = smoothInt*ETparams.samplingFreq;
assert(mod(smoothInt*ETparams.samplingFreq,1)==0,'window duration %.1f ms is not multiple of sample interval %.1f ms',smoothInt*1000,1000/ETparams.samplingFreq)

% Pixel values, velocities, and accelerations
%--------------------------------------------------------------------------
N = 2;                  % Order of polynomial fit
F = span;               % Window length (in samples)
assert(mod(span,2)==1,'Window duration %.1f ms should correspond to an uneven number of samples, not %d samples',smoothInt*1000, span);
fprintf('Window duration used: %1$.1f ms, %2$d samples (%3$d+1+%3$d)\n',smoothInt*1000,span,(span-1)/2);
[~,g] = sgolay(N,F);   % Calculate S-G coefficients

% Extract relevant gaze coordinates for the current trial.
X = data.X;
Y = data.Y;

% Calculate the velocity and acceleration
data.X = filter(g(:,1),1,X);
data.Y = filter(g(:,1),1,Y);

data.velX = filter(g(:,2),1,X);
data.velY = filter(g(:,2),1,Y);
data.vel = hypot(data.velX, data.velY)*ETparams.samplingFreq;

data.accX = filter(g(:,3),1,X);
data.accY = filter(g(:,3),1,Y);
data.acc = hypot(data.accX, data.accY)*ETparams.samplingFreq^2;