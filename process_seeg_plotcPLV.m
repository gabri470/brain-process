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
					distMask		= createSameReferenceMask(channels);
					
					[PLV, kPLV,iPLV,kiPLV] = computeAverageSingleSubject(data,channels,distMask);

					groupPLV(subjIdx,fIdx)  = PLV;
					groupKPLV(subjIdx,fIdx) = kPLV;
					groupiPLV(subjIdx,fIdx) = iPLV;
					groupKiPLV(subjIdx,fIdx)= kiPLV;

			end	
	end
	
	figure,
	subplot(2,2,1)
	semilogx(f,mean(groupPLV,1));
	xlim([min(f), max(f)]);
	ylabel('PLV')
	xlabel('Frequency');
	subplot(2,2,2)
	semilogx(f,mean(groupKPLV,1));
	ylabel('K')
	xlabel('Frequency');
	xlim([min(f), max(f)]);
	subplot(2,2,3)
	semilogx(f,mean(groupiPLV,1));
	xlim([min(f), max(f)]);
	ylabel('iPLV')
	xlabel('Frequency');
	subplot(2,2,4)
	semilogx(f,mean(groupKiPLV,1));
	ylabel('K')
	xlabel('Frequency');
	xlim([min(f), max(f)]);


	OutputFiles = {};

end

function [PLV,iPLV,kPLV,kiPLV] = computeAverageSingleSubject(data,channels,mask)
%COMPUTEAVERAGESINGLESUBJECT Description
%	ADJ = COMPUTEAVERAGESINGLESUBJECT(DATAIN) Long description
%
			nChans 			= numel(channels.Channel);

			cPLV				= complex(zeros(nChans),zeros(nChans));

			mask				= reshape(1:(nChans^2),[nChans,nChans]);

			cPLV(triu(mask)~=0) = (squeeze(data.TF(:,:,1)));

			cPLV(diag(mask)) = 0;

			PLVSurrMean = mean(squeeze(abs(data.TF(:,:,2))));

			iPLVConfLimit = std(squeeze(abs(imag(data.TF(:,:,2)))))*2.5758;
			PLVConfLimit = PLVSurrMean*2.42;

			PLV = abs(cPLV);
			iPLV= abs(imag(cPLV));

			PLV(PLV <= PLVConfLimit) = 0;
			iPLV(iPLV <= iPLVConfLimit) = 0;

			PLV(~mask) = 0;
			iPLV(~mask) = 0;


			kPLV = sum(sum(((PLV > PLVConfLimit))))/nchoosek(nChans,2);
			kiPLV= sum(sum(((iPLV > iPLVConfLimit)))) / nchoosek(nChans,2);

		  PLV = mean((PLV( PLV ~= 0 )));
			iPLV= mean((iPLV( iPLV ~= 0)));

end

function mask = createSameReferenceMask(channels)
%CREATESAMEREFERENCEMASK Description
%	MASK = CREATESAMEREFERENCEMASK(CHANNELS) Long description
%

channelNames = [{channels.Channel.Name}];

nChans = numel(channelNames);
mask = ones(nChans);

for channelIdx = 1:nChans
	currLabels = strsplit(channelNames{channelIdx},'-');

	rowMask = ~cellfun(@isempty,regexp(channelNames,regexprep(...
			(strjoin(['(^', currLabels(1),'-|-',currLabels(1),'$)'])),'\s','')));
	mask(channelIdx,rowMask) = 0;

	rowMask = ~cellfun(@isempty,regexp(channelNames,regexprep(...
			(strjoin(['(^', currLabels(2),'-|-',currLabels(2),'$)'])),'\s','')));

	mask(channelIdx,rowMask) = 0;

end
	

end
