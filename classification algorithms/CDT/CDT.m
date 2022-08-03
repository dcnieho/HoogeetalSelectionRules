% This is a stripped and modified version of the CDT algorithm, as used for
% Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. & Hessels,
% R.S. (2022). Fixation classification: how to merge and select fixation
% candidates. Behavior Research Methods.
% https://doi.org/10.3758/s13428-021-01723-1
%
% If you use this code, please cite the original paper Veneri et al. (2011)
% and the Hooge et al. (2022) paper for which the modified version was
% developed:
%
% Veneri, G., Piu, P., Rosini, F., Federighi, P., Federico, A., & Rufa, A.
% (2011). Automatic eye fixations identification based on analysis of
% variance and covariance. Pattern Recognition Letters, 32, 1588–1593.
% https://doi.org/10.1016/j.patrec.2011.06.012
%
% and
% 
% Hooge, I.T.C., Niehorster, D.C., Nyström, M., Andersson, R. & Hessels,
% R.S. (2022). Fixation classification: how to merge and select
% fixation candidates. Behavior Research Methods.
% https://doi.org/10.3758/s13428-021-01723-1

% USAGE: [list thoptim xi_w fv rho_w]=extract_fixation_ANOVA(x,y, window, th)
% where x and y are vector of the same length
% window = 6 (thats 25ms/acquisition rate) and th=0.05
%
% Based on F-DT
% Giacomo Veneri and Pietro Piu and Pamela Federighi and Francesca Rosini and Antonio Federico and Alessandra Rufa,
% "Eye Fixations Identication based on Statistical Analysis - Case study",
% 2010 IAPR Workshop on Cognitive Information Processing
%
% Modified on C-DT
% Giacomo Veneri and Pietro Piu and Francesca Rosini and Pamela Federighi  and Antonio Federico and Alessandra Rufa,
% "Automatic Eye Fixations Identification based on Analysis of Variance and Covariance", pattern Recognition Letters, 2011
% Please cite: Veneri & Piu 2010
%
% The algorithm is based on this assumption, and tests the hypothesis by a
% statistical method such as the F-test; the F-test is used to verify that
% two populations (with normal distribution) have the same variance (this is
% the so called "null hypothesis" H0 to be tested against the alternative
% complementary hypothesis H1 that the two populations have heterogeneous
% variances) and is a standard statistical procedure.
%
% return 
% list = [start fixation, end fixation, centroid of fixation, h ,p of f-test ]
% thoptim estimated threshold


% function [list thoptim xi_w fv rho_w]=veneri(x,y,params)
function fix=CDT(dat,channel,parameters)


% params
window      = round(parameters.windowLength/1000*dat.freq);
mergeWindow = 0;
thresh      = parameters.threshold;

list=[];
x = dat.(channel).X;
y = dat.(channel).Y;
%avoid mistake
nd=min(length(x),length(y));

% evaluate covariance
[xi_w, thoptim, fv, rho_w]=covWindow(x,y,window,thresh);
% sprintf('th: %f from %f',thoptim,thresh)

thoptim=max(thoptim,0.001);
%% look for centroid
pctv=1;

ppIdx = find(xi_w<thoptim);
pIdx = erode(ppIdx,window);

[npIdx, mpIdx]=size(pIdx);
if (npIdx>0)
%     sprintf('Found %d centroid',npIdx)
else
    sprintf('Found 0 centroid')
    show(x,y,xi_w,window,thoptim);
    return;
end

%% look for margin
list=[];
for ifix=1:npIdx,
    leftIdx = pIdx(ifix,2);
    sWin=max(1,leftIdx)-pIdx(ifix);
    rightIdx=pIdx(ifix,3);
    eWin=min(nd, rightIdx)-pIdx(ifix);
    list=[list ; [pIdx(ifix), sWin, eWin]];
end;


%% look for border
[n m] = size(list);
for ifix=1:n,
    winExtLeft=round(window/2);
    winExtRight=round(window/2);
    while((winExtRight>1) || (winExtLeft>1)),
        fixStart = max(1,list(ifix,1) + list(ifix,2));
        fixEnd = min(nd,list(ifix,1) + list(ifix,3));

        rv=nanvar(x(fixStart:fixEnd))*nanvar(y(fixStart:fixEnd));
        
        %look for intersection
        fixationsIdx=[list(:,1)+list(:,2) , list(:,1)+list(:,3)];

        % avoid intersection
        if (winExtRight>0) && (fixEnd+winExtRight<=nd) && isempty(find((fixationsIdx(:,1)<fixEnd+winExtRight) & (fixationsIdx(:,2)>fixEnd+winExtRight))),
            
            rvl=nanvar(x(fixStart:fixEnd+winExtRight)) * nanvar(y(fixStart:fixEnd+winExtRight));
            
            if (rvl<=rv*pctv)
                list(ifix,3) = list(ifix,3)+winExtRight;
            else
                winExtRight=winExtRight-1;
            end
        else
            winExtRight=0;
        end;

        % avoid intersection
        if (winExtLeft>0) && (fixStart-winExtLeft>0) && isempty(find((fixationsIdx(:,1)<fixStart-winExtLeft) & (fixationsIdx(:,2)>fixStart-winExtLeft))),

             rvl=nanvar(x(fixStart-winExtLeft:fixEnd)) * nanvar(y(fixStart-winExtLeft:fixEnd));
            
            if (rvl<rv*pctv)
                list(ifix,2) = list(ifix,2)-winExtLeft;
            else
                winExtLeft=winExtLeft-1;
            end
        else
            winExtLeft=0;
        end;
    end;
end;
%% adjust bounding box
list=round([max(list(:,1)+list(:,2),1) , min(list(:,1)+list(:,3),length(x)), list(:,1)]);


%% avoid over segementation
[n m]=size(list);
listMerged=[];
previousIdxEnd=-1;

for i=1:n,
    if (i>1) && (list(i,1) - previousIdxEnd<=mergeWindow) % merging
        [nm mm]=size(listMerged);
        listMerged(nm,2)=list(i,2); % update end
        listMerged(nm,3)=round((listMerged(nm,1)+listMerged(nm,2))/2); %update centroid
        [h,p] =vartest2(x(list(i,1):list(i,2)),y(list(i,1):list(i,2)));
        if isnan(h)
            h=0;
            p=1;
        end
        listMerged(nm,4:5) = [h,p];
        previousIdxEnd=listMerged(nm,2);
    else
        [h,p] =vartest2(x(list(i,1):list(i,2)),y(list(i,1):list(i,2)));
        if isnan(h)
            h=0;
            p=1;
        end
        listMerged=[listMerged; [list(i,:), [h, p]]];
        previousIdxEnd=list(i,2);
    end
end

%% return value <<bye bye>>
list=listMerged(listMerged(:,2)-listMerged(:,1)>1,:);

% % create the 'events' vector
% events = zeros(length(x),1);
% for i = 1:size(list,1)
%     events(list(i,1):list(i,2)) = 1; % Add 1s from the onset to the offset of each fixation
% end 

% prepare output
fix.start  = list(:,1);
fix.end    = list(:,2);

% this algo may deliver events with fixation starts or ends during missing
% data
% replace start in missing with end of first missing during that fixation
% replace end in missing with start of last missing during that fixation
[mon,moff] = bool2bounds(isnan(x));

while true
    qDoneSomething = false;
    for p=length(fix.start):-1:1
        qGat = mon>fix.start(p) & mon<fix.end(p) | moff>fix.start(p) & moff<fix.end(p);
        if any(qGat)
            % there is one or more holes in this fixation
            % we deal with one hole at a time. If there is more than one,
            % it is picked up on the next iteration
            iGat = find(qGat,1);
            hon  =  mon(iGat);
            hoff = moff(iGat);
            
            if hon>fix.start(p) && hoff<fix.end(p)
                % if hole(s) inside fixation -> split fixation into two
                fix.start = [fix.start(1:p); hoff+1; fix.start(p+1:end)];
                fix.end   = [fix.end(1:p-1);  hon-1; fix.end(p:end)];
            elseif hon>fix.start(p)
                % fixation starts in data, ends in hole -> move end to start of
                % hole
                fix.end(p) = hon-1;
            elseif hoff<fix.end(p)
                % fixation starts in hole, ends in data -> move start to end of
                % hole
                fix.start(p) = hoff+1;
            else
                % fixation starts in hole and ends in hole, so there's
                % no data during fixation -> remove
                fix.start(p) = [];
                fix.end(p) = [];
            end
            
            qDoneSomething = true;
        end
    end
    if ~qDoneSomething
        break;
    end
end

end

function void=show(x,y,cvn,w,th)
h=figure();
subplot(2,1,1);grid on;
hold on;
plot(x,'g');
plot(y,'b');
subplot(2,1,2);grid on;
hold on;
plot(cvn,'k');
plot(ones(length(cvn),1)*th,'r');
beep();
pause(5);
close(h);

end

function v=erode(x,window)
v=[];
if (length(x)<=0)
    return;
end
x=[-10*window; x];
dx = diff(x);
interruptedX = find(dx>1);

if (isempty(interruptedX)),
    v=[round(length(x)/2) min(x) max(x)];
end;

interruptedX=interruptedX+1;
for i=1:length(interruptedX)
    if (i+1>length(interruptedX)),
        next=length(x);
    else
        next = interruptedX(i+1)-1;
    end;
    c_min = min(x(interruptedX(i):next));
    c = mean(x(interruptedX(i):next));
    c_max = max(x(interruptedX(i):next));
    v = [v; c c_min c_max];
end;

v=round(v);
end

function v=varWindow(x,w)

[n m] = size(x);
v=[];

for i=1:n,
    is=max(1,i-w);
    in=min(n,i+w);
    v(i)=var(x(is:in));
end;
end

function v=fWindow(x,y,w)
[n m] = size(x);
v=zeros(n,1);

for i=1:n,
    is=max(1,i-w);
    in=min(n,i+w);
    [H,P] =vartest2(x(is:in),y(is:in));
    v(i)= P;
end;
end

function v=ncov(x,y)
v=nancov(x,y);
v=abs(v(1,2));
end

function [v th fv rho_w]=covWindow(x,y,w,thc)

[n m] = size(x);
v=zeros(n,1);

for i=1:n,
    is=max(1,i-w);
    in=min(n,i+w);
    v(i)=ncov(x(is:in),y(is:in));
end;

v=v./sqrt(nanvar(x)*nanvar(y));
rho_w=v;

%look for sure fixation => covariance is min
fv=fWindow(x,y,w);
idx=find(fv>0.05 & v<thc);
v(idx)=0;

%amplify small covariance
idx=find(fv<0.01 & v>=thc);
v(idx)=1;

% adjust threshold
idx=find(fv>0.05);
th=mean(v(idx));
if isnan(th)
    th=thc;
    sprintf('WARNING: too few samples')
end
end


