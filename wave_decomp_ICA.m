function [SWA, SWD, THR NSAMPLES] = wave_decomp_ICA(comp,N,wfilter)

[nchans, nsamples] = size(comp);

L = ceil(nsamples/(2^N))*(2^N);

comp(:,nsamples+1:L) = zeros(nchans,L-nsamples);

for j=1:size(comp,1)

    [thr,sohr,keepapp] = ddencmp('den','wv',comp(j,:));
    
    [swa,swd] = swt(comp(j,:),N,wfilter);
    
    %         SWD(abs(SWD)<THR) = 0;
    
    %         comp.trial{i}(j,:) = comp.trial{i}(j,:) - iswt(SWA,SWD,wfilter);
    SWD(:,:,j) = swd(:,1:nsamples);
    SWA(:,:,j) = swa(:,1:nsamples);
    THR(j) = thr;
    
end

    NSAMPLES = nsamples;
