function varargout = process_seeg_testCohort( varargin )

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
varargout =  {};
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Phase Synch';
sProcess.FileTag     = '__';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'SEEG';
sProcess.Index       = 901;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data'};
sProcess.OutputTypes = {'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(~, sInputs) %#ok<DEFNU>

	% init output PLV
	DataMat = db_template('timefreqmat');


	for fileIdx = 1:numel(sInputs) 
			% read one file at time
			data = in_bst_data(sInputs(fileIdx).FileName);
			channels = in_bst_channel(sInputs(fileIdx).ChannelFile);

			% [plv, nPLV, amp] = computePhaseMetric(data,nPermutation) 
			edgeData = computePhaseMetric(data, 3,'plv');

			% ===== SAVE THE RESULTS =====   
			% % Get the output study (pick the one from the first file)
			DataMat.TF				 	= edgeData.plv;

			iStudy 							= sInputs(fileIdx).iStudy;    
			DataMat.Comment     = 'cPLV';    
			DataMat.ChannelFlag = data.ChannelFlag;% List of good/bad channels (1=good, -1=bad)    
			DataMat.Time				= 1;
			DataMat.Freqs				= [1 2];
			DataMat.DataType    = 'data';    
			DataMat.Method      = 'plv';    
			DataMat.DataFile    = sInputs(fileIdx).FileName;    
			DataMat.Measure     = 'other';   
			DataMat.RefRowNames = [{channels.Channel.Name}];    
			DataMat.RowNames    = [{channels.Channel.Name}];    
			DataMat.nAvg        = length(sInputs); % Number of epochs that were averaged to get this file    
			% Create a default output filename     
			OutputFiles{1} = bst_process('GetNewFilename', ...
				fileparts(sInputs(fileIdx).FileName), 'timefreq_connectn_plv');    

			% Save on disk    
			save(OutputFiles{1}, '-struct', 'DataMat');    

			% Register in database    
			db_add_data(iStudy, OutputFiles{1}, DataMat);
	end

end

function [edgeData] = computePhaseMetric(data,nPermutation,metric)
% Description
%	[PTE,PLV, F] = computePhaseMetric(data,channelNames)

switch metric
		case 'plv'
				edgeData.plv = computePlv(data.F,nPermutation);
		case 'wpli'
				edgeData.wpli = [];
		case 'pTE'
				edgeData.pTE = [];
		otherwise
				bst_error(sprintf('metric %s unsupported',metric'));
end

end

function [out] = computePlv(signals,nPermutation)
% Description
%	[PLV,PLVSURR,AMPCORR,AMPCORRSURR] = computePlv(signals,nPermutation)

% Inst. phase [ch x T]
Xfc = hilbert(signals')';
% normalize 
Xfc = Xfc./abs(Xfc);

[nChans,nSamples] = size(Xfc);

% allocate fullmat plv mat 
cPlv = complex(eye(nChans),zeros(nChans));
cPlvSurr = complex(eye(nChans),zeros(nChans));

%compute only lower triangular 
for iCh = 1:nChans-1
		for jCh = iCh+1:nChans
				linIdx = sub2ind([nChans,nChans],iCh,jCh);
				cPlv(linIdx) = (Xfc(iCh,:)*Xfc(jCh,:)')./nSamples;
				offset = randi(nSamples,1);
				cPlvSurr(linIdx) = (Xfc(iCh,:)*circshift(Xfc(jCh,:),offset)')./nSamples;
		end
end


%% FASTER CODE? %%
%parfor linIdx = 2:nEdges+1
%	[iCh, jCh] = ind2sub([nChans,nChans],linIdx);
%	cPlv(linIdx) = (Xfc(iCh,:)*Xfc(jCh,:)')./nSamples;
%  offset = randi(nSamples,1);
%	cPlvSurr(linIdx) = (Xfc(iCh,:)*circshift(Xfc(jCh,:),offset)')./nSamples;
%end

% get vectorized mat
out(:,1,1) = cPlv(cPlv ~= 0);
out(:,1,2) = cPlvSurr(cPlv ~= 0); 

end
