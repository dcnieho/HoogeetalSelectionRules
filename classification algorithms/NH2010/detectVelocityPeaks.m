function data = detectVelocityPeaks(data)
 

 % Find peaks larger than a threshold ('ETparams.peakDetectionThreshold')
 % Sets a '1' where the velocity is larger than the threshold and '0'
 % otherwise
data.InitialVelPeakIdx  = (data.vel > data.peakDetectionThreshold);
 
 
