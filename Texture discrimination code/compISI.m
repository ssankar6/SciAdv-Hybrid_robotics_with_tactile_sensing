function [ISIm,ISI]=compISI(v,buffer)
% % compute ISI for a given window
% inputs
% % v is an array of membrane potential values where the rows are different neurons and the columns are time stamps of the neuron's membrane potential
% % buffer is the max amount of time to wait (in ms) for another spike

% outputs
% % ISI is a cell array containing all the ISIs generated for each sensor
% % ISIm is an array of average ISIs since the ISI distribution will vary between neurons

for i=1:size(v,1)
    spikeidx=find(v(i,:)==30); % generate time stamps of spikes
    delta=diff(spikeidx); % difference in time stamps -> ISI
    for j=1:length(delta)
        if delta(j)>buffer % longer than a second...
            delta(j)=NaN; % ISI not valid since no apparent stimulus in dead time
        end
    end
    delta = delta';
    delta = delta(~isnan(delta))'; % find actual ISIs
    ISI{i}=delta;
    if isempty(ISI{i})
        ISI{i}=0; % no spikes -> no ISI
    end
    if isnan(mean(delta))
        ISIm(i)=0;
    else
        ISIm(i)=mean(delta); % ms
    end
end
end