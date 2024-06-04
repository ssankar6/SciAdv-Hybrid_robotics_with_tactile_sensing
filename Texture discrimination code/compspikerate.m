function sr = compspikerate(v,bin)
% % compute spike rates in bins
% inputs
% % v is an array of membrane potential values where the rows are different neurons and the columns are time stamps of the neuron's membrane potential
% % bin is the length of bins (in ms) over a time window to analyze data in

% output - sr is an array of spike rates

numbins= floor(size(v,2)/bin); % number of bins based on bin length and duration of overall neuron dynamics
sr=zeros(size(v,1),numbins);
for i=1:numbins
    spikecount=v(:,[i-1]*bin+1:[i-1]*bin+bin)>=30; % find how many spikes occur
    sr(:,i)=sum(spikecount,2)/bin;
end
end