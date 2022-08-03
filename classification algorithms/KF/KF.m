% This is a stripped and modified version of the KF algorithm, as used
% for Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. &
% Hessels, R.S. (2022). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.
% https://doi.org/10.3758/s13428-021-01723-1
%
% If you use this code, please cite the original paper Komogortsev et al.
% (2010) and the Hooge et al. (2022) paper for which the modified
% version was developed:
%
% Komogortsev, O. V., Gobert, D. V., Jayarathna, S., Koh, D. H., & Gowda,
% S. M. (2010). Standardization of automated analyses of oculomotor
% fixation and saccadic behaviors. IEEE Transactions on Biomedical
% Engineering, 57(11), 2635–2645.
%
% and
% 
% Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. & Hessels,
% R.S. (2022). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.
% https://doi.org/10.3758/s13428-021-01723-1

function evts = KF(dat,channel,parameters,which_event)
KFchi_threshold = parameters.chi_threshold;
KFwindow_size   = parameters.window_size;
KFdeviation     = parameters.deviation;

nSamp       = length(dat.(channel).X);
delta_t_sec = 1/dat.freq;
qMiss       = isnan(dat.(channel).X);
eye_records = nan(nSamp,6);
eye_records(:,1) = dat.(channel).X(:);
eye_records(:,2) = dat.(channel).Y(:);
% eye_records(:,3) = nan;
eye_records(qMiss,4) = 4; %obj.NOISE_TYPE;
eye_records(:,5) = double(~qMiss);
eye_records(:,6) = 1:length(dat.(channel).X);

% Kalman filter setup
KF_K = 0;                          % This is Kalman filter's variables
KF_x = [0; 0];                     % For any details look for Kalman filter information
KF_y = [0; 0];                     % ...
KF_P = [1 0; 0 1;];                % ...
KF_A = [1 delta_t_sec; 0 1];       % ...
KF_H = [1 0];                      % ...
KF_Q = [0.5 0; 0 0.5];             % ...
KF_R = 0.5;                        % ...

KF_predicted = zeros(nSamp,1); % This array holds the predicted speeds
KF_measured = zeros(nSamp,1);  % This array holds the predicted speeds
KF_x_velocity_degree = zeros(nSamp,1);
KF_y_velocity_degree = zeros(nSamp,1);
% Calculate absolute degree velocity of our records
KF_x_velocity_degree( 2:end ) =(    eye_records( 2:end,1 ) - ...
    eye_records( 1:end-1,1) ) / delta_t_sec;
KF_y_velocity_degree( 2:end ) =(    eye_records( 2:end,2 ) - ...
    eye_records( 1:end-1,2 ) ) / delta_t_sec;
% First point is a special case
KF_x_velocity_degree(1) = 0;
KF_y_velocity_degree(1) = 0;
eye_records(:,3) = sqrt( KF_x_velocity_degree.^2 + KF_y_velocity_degree.^2 );
% First point is a special case
eye_records(1,3) = 0;

% Generating predicteded velocities using Kalman filter
for i=1:length(eye_records)
    % Discarding the noise
    if ( eye_records(i,5) == 0 || isnan(eye_records(i,3)) )
        eye_records(i,4) = 4; %obj.NOISE_TYPE;
    else
        % Prediction step of Kalman filter
        KF_x = KF_A * KF_x;
        KF_y = KF_A * KF_y;
        KF_P = KF_A * KF_P * KF_A.' + KF_Q;
        % Calculate predicteded velocity, storing current predicteded and measured velocities
        KF_predicted(i) =   sqrt(KF_x(2).^2 + KF_y(2).^2);
        KF_measured(i) =    eye_records(i,3);
        % Update step of Kalman filter
        % Cory code
        KF_K = (KF_P * KF_H.') / (KF_H * KF_P * KF_H.' + KF_R);
        % Sam's code
        %                    KF_K = KF_P.' * KF_H.' / (KF_H * KF_P.' * KF_H.' + KF_R);
        KF_x = KF_x + KF_K * (KF_H * [eye_records(i,1); KF_x_velocity_degree(i)] - KF_H * KF_x);
        KF_y = KF_y + KF_K * (KF_H * [eye_records(i,2); KF_y_velocity_degree(i)] - KF_H * KF_y);
        KF_P = KF_P - KF_K * KF_H * KF_P;
    end
end

for i=1:length(eye_records)
    KF_window_start = i;
    KF_window_ended = i+KFwindow_size;
    if(KF_window_ended > length(eye_records) )
        break;
    end
    KF_chi=0;
    for j=KF_window_start:KF_window_ended
        if( eye_records( j,5) == 1 )
            KF_chi = KF_chi + ( KF_measured(j) - KF_predicted(j) ).^2 / KFdeviation;
        end
    end
    if( eye_records( KF_window_ended,5 ) == 1)
        if (abs(KF_chi) < KFchi_threshold)
            % Classify points below Chi-square threshold as fixations
            eye_records(KF_window_ended,4) = 1; %FIXATION_TYPE;
        else
            % Classify points above Chi-square threshold as saccades
            eye_records(KF_window_ended,4) = 2; %SACCADE_TYPE;
            % Due to our valocity calculations (left difference) we should include
            % previous point to the saccade too
            if( KF_window_ended > 1)
                if( eye_records( KF_window_ended - 1,5 ) == 1 )
                    eye_records( KF_window_ended - 1,4 ) = 2; %SACCADE_TYPE;
                end
            end
        end
    end
end

% Special case (unclassified -> noise)
eye_records( ( isnan(eye_records(:,4)) ),4 ) = 4; %NOISE_TYPE;



% prepare output
assert(ismember(which_event,[1 2])) % should be fix (1) or sac (2)
[evts.start,evts.end] = bool2bounds(eye_records(:,4)==which_event);
