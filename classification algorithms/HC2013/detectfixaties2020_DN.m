function [fmark,thr2,meanvel,stdvel] = detectfixaties2020_DN(mvel,f)

% opschoonactie
% 16 oktober 2011 IH
% 14 januari 2020 IH
% 01 march 2021 DN: keep only core algorithm

thr             = f.thr;
tc              = f.tc;

qvel            = mvel < thr;                      % kijk waar de snelheid kleiner is dan thr
meanvel         = mean(mvel(qvel));                % bepaal mean van de velocity tijdens de fixaties
stdvel          = std(mvel(qvel));                 % bepaal de std van de velocity tijdens de fixaties

counter         = 0;
oldthr          = 0;
while 1
    thr2        = meanvel + f.lambda*stdvel;
    qvel        = mvel < thr2;                     % kijk waar de snelheid kleiner is dan thr
    
    if round(thr2) == round(oldthr) || counter == f.counter % f.counter voor maximaal aantal iteraties
        break;
    end
    meanvel     = mean(mvel(qvel));
    stdvel      = std(mvel(qvel));                 % bepaal de std van de velocity tijdens de fixaties    
    oldthr      = thr2;
    counter     = counter + 1;
end

thr2            = meanvel + f.lambda*stdvel;       % bepaal nieuwe drempel gebaseerd op de ruis in de data
qvel            = mvel < thr2;                     % kijk waar de snelheid kleiner is dan thr
[on,off]        = detectswitches(qvel');           % bepaal fixaties

% we willen fixatie einden en -startpunten bepalen. Die -1 en +1 beneden
% lijken verwarrend, maar als tc (sample interval) groot is (groot t.o.v.
% saccadeduur) kun je geen accurate saccadestart en -eindpunten bepalen en
% wil je zo min mogelijk fixatie weggooien.

on(2:end)       = on(2:end) - 1;                   % on is de fixatiestart (index, dus in samples, niet in tijd). Met de eerste doe ik niets want die staat al goed
off(1:end-1)    = off(1:end-1) + 1;                % off is het einde van de fixatie (index, dus in samples, niet in tijd). Met de laatste doe ik ook niet want ie staat ook goed
fmark           = sort([on off]);                  % gooi de markers op volgorde achter elkaar
