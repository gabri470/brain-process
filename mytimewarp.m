% timewarp()   - Given two event marker vectors, computes a matrix
%                that can be used to warp a time series so that its
%                evlatencies match newlatencies. Values of the warped
%                timeserie that falls between two frames in the original
%                timeserie will be linear interpolated.
% Usage:
%   >> warpmat = timewarp(evlatency, newlatency)
%
% Necessary inputs:
%   evlatency  - [vector] event markers in the original time series, in frames
%                Markers must be ordered by increasing frame latency.
%                If you want to warp the entire time-series, make sure
%                the first (1) and last frames are in the vector.
%   newlatency - [vector] desired warped event time latencies. The original
%                time series will be warped so that the frame numbers of its
%                events (see evlatency above) match the frame numbers in
%                newlatency. newlatency frames must be sorted by ascending
%                latency. Both vectors must be the same length.
%
% Optional outputs:
%      warpmat - [matrix] Multiplying this matrix with the original
%                time series (column) vectors performs the warping.
%
% Example:
%      % In 10-frame vectors, warp frames 3 and 5 to frames 4 and 8,
%      % respectively. Generate a matrix to warp data epochs in
%      % variable 'data' of size (10,k)
%      >> warpmat = timewarp([1 3 5  10], [1 4 8 10])
%      >> warped_data = warpmat*data;
%
% Authors: Jean Hausser, SCCN/INC/UCSD, 2006
%
% See also: angtimewarp(), phasecoher(), erpimage()
%
function M=mytimewarp(newLatency, oldLatency,ord)

if min(sort(oldLatency) == oldLatency) == 0
    error('oldLatency should be in ascending order');
    return;
end
if min(sort(newLatency) == newLatency) == 0
    error('newLatency should be in ascending order');
    return;
end
if length(oldLatency) ~= length(newLatency)
    error('oldLatency and newLatency must have the same length.');
    return;
end
if length(oldLatency) < 2 | length(newLatency) < 2
    error(['There should be at least two events in evlatency and ' ...
        'newlatency (e.g., "begin" and "end")'] );
    return;
end
if oldLatency(1) ~= 1
    disp(['Assuming old and new time series beginnings are synchronized.']);
    disp(['Make sure you have defined an ending event in both the old and new time series!']);
    oldLatency(end+1)=1;
    newLatency(end+1)=1;
    oldLatency = sort(oldLatency);
    newLatency = sort(newLatency);
end

t = 1:max(oldLatency);

x = piecewise_linear_func(t,oldLatency,newLatency);
% x = lagrange_pol_evaluation(t,newLatency,oldLatency,1)

% sze = numel(x);
% 
% xa = floor(x);
% xb = floor(x)+1;
% wb = (x-xa);
% wa = (xb-x);
% 
% M = (zeros(sze,sze));
% 
% M([1:sze]+(xa-1)*sze) = wa;
% M([1:sze-1]+(xb(1:end-1)-1)*sze) = wb(1:end-1);
% 
M = lagrange_pol_matrix(x,ord);

end

function y = piecewise_linear_func(t,To_latency,From_latency)

sze = numel(t);
num_evt = numel(From_latency);

y = zeros(1,sze);

for i=1:num_evt-1
    y(From_latency(i):From_latency(i+1)) = (t(From_latency(i):From_latency(i+1))-From_latency(i)).*(To_latency(i+1)-To_latency(i))./...
        (From_latency(i+1)-From_latency(i))+...
        To_latency(i);
end

end

function val = lagrange_pol_matrix(x,ord)

sze = numel(x);
x = [x sze+[1:ord-1]]';
floor_x = bsxfun(@plus,floor(x),[0:(ord-1)]);
val = zeros(sze,sze+ord-1);
for i = 1:ord
    val([1:sze]'+(floor_x(1:sze,i)-1)*sze) = lagrange_coeff(x(1:sze),floor_x(1:sze,setdiff(1:ord,i)),floor_x(1:sze,i));
end
%     val = lagrange_pol_matrix_coef(x,bsxfun(@plus,floor(x),[0:(ord-1)]),1);

val(:,end-(ord-1)+1:end) = [];
end

function val = lagrange_pol_matrix_coef(x,xk,j)

if(j==size(xk,2)+1)
    val = [];
else
    val = [lagrange_coeff(x,xk(:,setdiff(1:size(xk,2),j)),xk(:,j)) lagrange_pol_matrix_coef(x,xk,j+1)];
end
end

function val = lagrange_pol_evaluation(x,xk,yk,j)

if(j==size(yk,2)+1)
    val = 0;
else
    val = yk(j)*lagrange_coeff(x,xk(:,setdiff(1:numel(xk),j)),xk(:,j))+lagrange_pol_evaluation(x,xk,yk,j+1);
end
end

function L = lagrange_coeff(x,xk,xj)

if(isempty(xk))
    L=1;
else
    L = (x-xk(:,1))./(xj-xk(:,1)).*lagrange_coeff(x,xk(:,2:end),xj);
end
end
