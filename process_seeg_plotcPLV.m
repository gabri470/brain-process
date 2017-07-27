function varargout = process_seeg_plotcPLV( varargin )

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
sProcess.Comment     = 'Plot cPlv afo F';
sProcess.FileTag     = '__';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'SEEG';
sProcess.Index       = 1001;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'timefreq'};
sProcess.OutputTypes = {'timefreq'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(~, sInputs) %#ok<DEFNU>
	
	subjectNames = [{sInputs.SubjectName}];
	nSubjects    = numel(unique(subjectNames));

  comments = cell(1,numel(sInputs));
	for fileIdx = 1:numel(sInputs)
				
			comments{fileIdx} = getfield(in_bst_data(sInputs(fileIdx).DataFile,'Comment'),'Comment');;

	end

	fStrings 			= regexp(comments,'\d+.*\d+','match');
	fStrings 			= [fStrings{:}];
	nominalFreqs 	= str2double(fStrings);
	f 			 			= unique(str2double(fStrings));
	nFreqs				= numel(f);

	for subjIdx = 1:nSubjects
			for fIdx = 1:nFreqs

			    subjectMask = ~cellfun(@isempty,regexp({sInputs.SubjectName},subjectNames{subjIdx}));
					fMask				= nominalFreqs == f(fIdx);
    			fileIdx 		= find(subjectMask & fMask);

					% read one file at time
					data 				= in_bst_timefreq(sInputs(fileIdx).FileName);
					channels 		= in_bst_channel(sInputs(fileIdx).ChannelFile);
					
					[adj, k] 		= computeAverageSingleSubject(data,channels);

					EdgeData(subjIdx,fIdx) = adj;
					K(subjIdx,fIdx) = k;
			end	
	end
	
	figure,
	subplot(2,1,1)
	semilogx(f,mean(EdgeData,1));
	ylabel('PLV')
	subplot(2,1,2)
	semilogx(f,mean(K,1));
	ylabel('K')
	xlabel('Frequency');



	OutputFiles{1} = 0;

end

function [adj,k] = computeAverageSingleSubject(data,channels)
%COMPUTEAVERAGESINGLESUBJECT Description
%	ADJ = COMPUTEAVERAGESINGLESUBJECT(DATAIN) Long description
%
			nChans 			= numel(channels.Channel);

%			cPLV				= complex(zeros(nChans),zeros(nChans));
%			cPLVSurr		= complex(zeros(nChans),zeros(nChans));
	
			cPLV				= zeros(nChans);

			mask				= reshape(1:(nChans^2),[nChans,nChans]);

			cPLV(triu(mask)~=0) = abs(squeeze(data.TF(:,:,1)));

			cPLV(diag(mask)) = 0;

			PLVSurrMean = mean(squeeze(abs(data.TF(:,:,2))));

			k = sum(sum(((abs(cPLV)./PLVSurrMean) > 2.42)))/nchoosek(nChans,2);

			cPLV((abs(cPLV)./PLVSurrMean) < 2.42) = 0;

		  adj = mean(abs(cPLV( cPLV ~= 0 )));

end
