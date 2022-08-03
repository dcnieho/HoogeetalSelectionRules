% This is a stripped and modified version of the MST algorithm, as used
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

function evts = MST(dat,channel,parameters,which_event)
saccade_detection_threshold = parameters.saccade_detection_threshold;
window_size                 = parameters.window_size;

% note that while code seems to suggest 200 samples was used, 200 ms makes
% more sense and is what is written in Komogortsev, Gobert et al 2010. The
% algo seems to perform better as well with this smaller window
window_size                 = round(window_size/1000*dat.freq);


nSamp       = length(dat.(channel).X);
qMiss       = isnan(dat.(channel).X);
eye_records = nan(nSamp,6);
eye_records(:,1) = dat.(channel).X(:);
eye_records(:,2) = dat.(channel).Y(:);
% eye_records(:,3) = nan;
eye_records(qMiss,4) = 4; %obj.NOISE_TYPE;
eye_records(:,5) = double(~qMiss);
eye_records(:,6) = 1:length(dat.(channel).X);

% Begin processing data records by moving windows
window_count = ceil(nSamp / window_size );
for k=1:window_count
    % Get begin, end, and length of the current window
    start_record = window_size * (k-1) + 1;
    ended_record = min(start_record + window_size - 1,nSamp);
    window_length = ended_record - start_record + 1;
    % Computing Euclidean distance between each eye position inside the window
    X = zeros(window_length,2);
    X(:,1) = eye_records(start_record:ended_record, 1);
    X(:,2) = eye_records(start_record:ended_record, 2);
    distance_matrix = squareform(pdist(X,'euclidean'));
    distance_matrix(logical(eye(size(distance_matrix)))) = inf;
    
    % Building minimum spanning tree using Prim's algorithm. I'm afraid that
    % this implementation is not acceptable for large data arrays. It's better
    % to use something to improve this.
    vertices = zeros( window_length,1 );    % List of visited vertices
    vertices(1) = 1;                        % Mark first vertex as visited
    count_v = 1;                            % Total amount of visited vertices
    while (count_v<window_length)           % While some of the vertices is unvisited
        min_distance = inf;                 % Searching for minimal distance betwenn visited and unvisited vertices
        [v1, v2] = deal([]);
        for i=1:window_length
            if (vertices(i) == 1)
                for j=1:window_length
                    if (vertices(j) == 0)
                        if( min_distance > distance_matrix( i,j ) )
                            min_distance = distance_matrix( i,j );
                            v1 = i;
                            v2 = j;
                        end
                    end
                end
            end
        end
        vertices(v2) = 1;                   % Mark selected unvisited vertex as visited
        % Comparing new edge of the MST with threshold value
        if( min_distance < saccade_detection_threshold )
            % If length of this edge is less than threshold value then fixation
            % detected but only for valid data
            if( eye_records( v1 + start_record -1,5) == 1)
                eye_records( v1 + start_record - 1 ,4 ) = 1; %FIXATION_TYPE;
            else
                eye_records( v1 + start_record - 1 ,4 ) = 4; %obj.NOISE_TYPE;
            end
            if( eye_records( v2 + start_record -1,5) == 1)
                eye_records( v2 + start_record - 1 ,4 ) = 1; %obj.FIXATION_TYPE;
            else
                eye_records( v2 + start_record - 1 ,4 ) = 4; %obj.NOISE_TYPE;
            end
        else
            % saccade detected but only for valid data
            if( eye_records( v1 + start_record -1,5) == 1)
                eye_records( v1 + start_record - 1 ,4 ) = 2; %obj.SACCADE_TYPE;
            else
                eye_records( v1 + start_record - 1 ,4 ) = 4; %obj.NOISE_TYPE;
            end
            if( eye_records( v2 + start_record -1,5) == 1)
                eye_records( v2 + start_record - 1 ,4 ) = 2; %obj.SACCADE_TYPE;
            else
                eye_records( v2 + start_record - 1 ,4 ) = 4; %obj.NOISE_TYPE;
            end
        end
        count_v = count_v + 1;
    end
end
% Due to our method of velocity computations (left difference) we have to
% add left edge to our saccades
tmp_type = eye_records(2:end,4);
tmp_type(length(tmp_type)+1) = NaN;
eye_records( ((eye_records(:,5) == 1) & (tmp_type(:) == 2)),4) = 2;

eye_records( isnan(eye_records( :, 4 )),4) = 4; % unclassified -> noise



% prepare output
assert(ismember(which_event,[1 2])) % should be fix (1) or sac (2)
[evts.start,evts.end] = bool2bounds(eye_records(:,4)==which_event);
