function An=fs_dict(A)
% % feature scale a dictionary to normalize it
% input - A is a dictionary where the rows are the trials and the columns are the features

% output - An is a normalized dictionary

sub=min(A,[],2);
div=max(A,[],2);
An=(A-sub)./div;
end