
function iPk = fast_findpeaks(Y,th)

% take the sign of the first sample derivative
s = [ones(size(Y,1),1) sign(diff(Y,[],2))];

% find local maxima
iMax = [(diff(s,[],2)<0) false(size(Y,1),1)];

iPk = (Y > th) & iMax;

end