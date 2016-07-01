function [fft_sig, f] = mypwelch(sig,wdw,nfft,fs)

f = linspace(0,fs/2,floor(nfft/2)+1);

if(numel(wdw)==1)
    wdw = hamming(wdw,'symmetric');
end

window_l = length(wdw);

tap_num = floor(size(sig,2)/window_l);

% remove DC component
sig = sig(:,1:tap_num*window_l);
sig = sig - repmat(mean(sig,2),[1,size(sig,2)]);

sig = reshape(sig,[size(sig,1),window_l,tap_num]);
U = wdw'*wdw;

wdw = repmat(wdw',[size(sig,1),1,tap_num]);

xx = fft(sig.*wdw,nfft,2)/window_l;

fft_sig = mean(xx.*conj(xx),3);      % Auto spectrum.

% for i = 1:tap_num
%     
%     x = sig((1:window_l)+(i-1).*window_l,:).*repmat(window,1,size(sig,2));
%     
%     xx = fft(x,nfft)/window_l;
%     
%     P = xx.*conj(xx);      % Auto spectrum.
%     
%     fft_sig = fft_sig + P;
% end
% 
% fft_sig = fft_sig/tap_num;

fft_sig = fft_sig(:,1:floor(nfft/2)+1);

fft_sig(:,2:end) = 4*fft_sig(:,2:end);
