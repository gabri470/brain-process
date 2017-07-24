function varargout = process_filter_fir( varargin )

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
sProcess.Comment     = 'filter w fir bank';
sProcess.FileTag     = @GetFileTag;
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'Pre-process';
sProcess.Index       = 61;

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','raw'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

sProcess.options.minFreq.Type	 					= 'text';
sProcess.options.minFreq.Comment 				= 'Min Freq';
sProcess.options.minFreq.Value	 				= '2';

sProcess.options.scalingFactor.Type	 		= 'text';
sProcess.options.scalingFactor.Comment 	= 'Scaling Factor';
sProcess.options.scalingFactor.Value	 	= num2str(sqrt(2));

sProcess.options.bandFlatTop.Type 			= 'text';
sProcess.options.bandFlatTop.Comment 		= 'Band Flat Top';
sProcess.options.bandFlatTop.Value 			= '2';

sProcess.options.stopBand.Type 					= 'text';
sProcess.options.stopBand.Comment 			= 'Stop Band';
sProcess.options.stopBand.Value 				= '2';

sProcess.options.bandStopWidth.Type 		= 'text';
sProcess.options.bandStopWidth.Comment 	= 'Band Stop Width';
sProcess.options.bandStopWidth.Value 		= '5';

sProcess.options.nFilters.Type 					= 'text';
sProcess.options.nFilters.Comment 			= 'Number of Filters';
sProcess.options.nFilters.Value 				= '13';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
	Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

		fs = 1000;
		[filterSettings,fn] = generateFilterSettings(sProcess.options,fs);
		nFilters = size(filterSettings,1);
		bst_progress('start','Filtering','Filtering');

		totElement = numel(sInputs)+nFilters;

		for fileIdx = 1:numel(sInputs)
				for fIdx = 1:nFilters

						bst_progress('inc',((fileIdx-1)*numel(sInputs))+fIdx);

						DataMat = in_bst_data(sInputs(fileIdx).FileName);

						fs = 1/mean(diff(DataMat.Time));

						[filteredDataMat] = Compute(DataMat.F,filterSettings,fIdx);

						% ===== SAVE THE RESULTS =====   
						% Get the output study (pick the one from the first file)
						iStudy 							= sInputs(fileIdx).iStudy;    
						DataMat.Comment     = strcat('band ',num2str(fn));    
						DataMat.F						= filteredDataMat;
						DataMat.DataType    = 'data';    
						DataMat.nAvg        = length(sInputs); % Number of epochs that were averaged to get this file    

						% Create a default output filename     
						OutputFiles{1} = bst_process('GetNewFilename', ...
							fileparts(sInputs(fileIdx).FileName), 'data_band');    
						% Save on disk    
						save(OutputFiles{1}, '-struct', 'DataMat');    

						% Register in database    
						db_add_data(iStudy, OutputFiles{1}, DataMat);
				end
		end

end


function data = Compute(data,filterSettings,fIdx)

		[nFilters,nHarmonics] = size(filterSettings);

		for hIdx = 1:nHarmonics
				if ~isempty(filterSettings{fIdx,hIdx})
								 
						filt = designfilt(filterSettings{fIdx,hIdx}{:});

						filtDelay = round(mean(grpdelay(filt)));
						
						data = [data, flipud(data(:,end-filtDelay+1:end))];
						data = (filter(filt, data'))';
						data = data(:,filtDelay+1:end);
				end

		end


end


function [filterSettings,f] = generateFilterSettings(options,fs) 
%GENERATEFILTERSETTINGS wrapper for designfilt 
%	FILTERSETTINGS = GENERATEFILTERSETTINGS(SPROCESS.OPTIONS,FS) 
%	This function should translate the sProcess.options into a cellarray of fields 
%	that can be easily passed to designfilt


	nFilters 	= str2num(options.nFilters.Value);
	fn 			 	= str2num(options.minFreq.Value);
	wb			 	= str2num(options.bandFlatTop.Value);
	ws			 	= str2num(options.stopBand.Value);
	m				 	= str2num(options.scalingFactor.Value);

	f 			 	= zeros(nFilters,1);
	fnyq		 	= fs/2;
	lineNoise	= 50;
	harmonics	= linspace(50,fs/2,(fs/2)/50)/(fs/2);
	tmp 			= linspace(50,fs/2,(fs/2)/50);
	nHarmonics= numel(harmonics);

	filterSettings = cell(nFilters,nHarmonics+1);

	for fIdx = 1:nFilters

			if fn > 100
					m = sqrt(sqrt(2));
			end

			fprintf('[%f - %f] -> %f <- [%f - %f]\t',(fn / ws),(fn - (wb * fn)/2),fn,...
					(fn + (wb * fn)/2),(min([0.95*(fs/2) fn * ws])));
	
			passBandLp = (fn + (wb * fn)/2)/(fs/2);
			passBandHp = (fn - (wb * fn)/2)/(fs/2);
			stopBandLp = (min([0.95*(fs/2) fn * ws]))/(fs/2);
			stopBandHp = (fn / ws)/(fs/2);

			f(fIdx) = fn;
			fn = fn * m;

			filterSettings{fIdx,end} = {'bandpassfir',...
					'StopbandFrequency1',stopBandHp,...
					'PassbandFrequency1',passBandHp,...
					'PassbandFrequency2',passBandLp,...
					'StopbandFrequency2',stopBandLp,...
					'StopbandAttenuation1',80,...
					'StopbandAttenuation2',80,...
					'DesignMethod','kaiserwin'};

			
			curHarmonics = find((harmonics > stopBandHp & harmonics < stopBandLp));
			harmonics = harmonics.*(fs/2);

			for harmIdx = curHarmonics
					fprintf(' %f ', tmp(harmIdx));
					filterSettings{fIdx,harmIdx} = {'bandstopfir',...
							'FilterOrder',1500, ...
		          'CutoffFrequency1',harmonics(harmIdx)-(str2num(options.bandStopWidth.Value)),...
							'CutoffFrequency2',harmonics(harmIdx)+(str2num(options.bandStopWidth.Value)), ...
							'SampleRate',fs};
					
			end
			fprintf('\n');

	end

end

