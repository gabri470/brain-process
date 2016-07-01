function artefact = compute_artefact(swd,swa,thr,nsample,wfilter)

swd(abs(swd)<thr) = 0;
swa(abs(swd)<thr) = 0;

N = size(swd,1);
L = ceil(nsample/(2^N))*(2^N);

swd(:,nsample+1:L) = zeros(N,L-nsample);
swa(:,nsample+1:L) = zeros(N,L-nsample);

artefact = iswt(swa,swd,wfilter);

artefact = artefact(1:nsample);