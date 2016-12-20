function varargout = process_PSDStatisticalSignficance( varargin )
% PROCESS_EXAMPLE_CUSTOMAVG: Example file that reads all the data files in input, and saves the average.

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
% Authors: Francois Tadel, 2013

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'PSD Statistical significance';
    sProcess.FileTag     = ' | significant';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Frequency';
    sProcess.Index       = 1000;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
    sProcess.options.pval.Comment = 'p-val';
    sProcess.options.pval.Type    = 'value';
    sProcess.options.pval.Value   = {5, 'units', 5};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Initialize returned list of files
    OutputFiles = {};
    % Get option values
    pval = sProcess.options.pval.Value{1};
    
    % Reading all the input files in a big matrix
		for i = 1:length(sInputs)

				DataMat = in_bst(sInputs(i).FileName, [], 0);
				parentStruct = bst_process('GetInputStruct',sInputs(i).DataFile); 
				chFlag = in_bst_data(parentStruct.FileName,'ChannelFlag');
				
				goodChannelMask = chFlag.ChannelFlag == 1;
				psds  = squeeze(DataMat.TF);
				H 	  = significance(psds,DataMat.Freqs,pval);
				H	  = H & goodChannelMask;

				fprintf('%s %d\n',sInputs(i).FileName,sum(H));
				% ===== SAVE THE RESULTS =====
				% Get the output study (pick the one from the first file)
				iStudy = sInputs(i).iStudy;
				% Create a new data file structure
				DataMat.TF			= permute(psds(H,:),[1 3 2]);
				DataMat.RowNames	= DataMat.RowNames(H);
				DataMat.Comment		= strcat(DataMat.Comment, sProcess.FileTag);
				DataMat.History		= [DataMat.History {datestr(now)} ...
					{'Statistical Significance Test'}];

				OutputFiles{1} = bst_process('GetNewFilename', ...
					fileparts(sInputs(i).FileName), 'timefreq_psd_signficance');

				% Save on disk
				save(OutputFiles{1}, '-struct', 'DataMat');
				% Register in database
				db_add_data(iStudy, OutputFiles{1}, DataMat);
		end
end

function H = significance(psds,f1,alpha)
%SIGNIFICANCE Description
%	[H] = SIGNIFICANCE(PSDS,F1,N,ALPHA) Long description
%

%		f  			= f1(f1>5 & f1<=90)';
%		psd_data 	= psds(:,f1>5 & f1<=90);
		f 			= f1;
		mask 		= (f >= 10 & f <= 30);

%		A 			= [log10(f)  ones(size(f))];
		% A*x = b
%		coeffs		= pinv(A)*log10(psd_data)';

		bound = [5 10 30 80];
		options = optimoptions(@fminunc,'Algorithm',...
			'trust-region','MaxIter',1500,...
			'TolFun',1e-12,'TolX',1e-12);

		% ch x f
%		noise = repmat(10.^coeffs(2,:),size(f,1),1).*...
%		(repmat(f,1,size(psds,1)).^repmat(coeffs(1,:),size(f,1),1));

		H = logical(zeros(size(psds,1),1));

		for ch = 1:size(psds,1)

				I = psds(ch,:);
				param = fminunc(@(param)functional(I,param,f,bound),...
					[3,10,3],options);

%				P = noise(:,ch)';
								
				P = (pink_spectrum(param,f));
				U = (2*(55) /  chi2inv(alpha,2*(55))).*P;

				I = I./sum(I);

%				H(ch) = sum(I(mask) >= U(mask)) >= 5;
				H(ch) = iterSumBool(I(mask) >= U(mask));

		end % channel loop

end

function count = iterSumBool(vectBool)
%ITERSUMBOOL Description
%	H = ITERSUMBOOL(VECTBOOL) Long description
%
	count = 0;

	for ii = 1:numel(vectBool)

		if vectBool(ii)
			count = count + 1;
		else
			count = 0;
		end
	
	end; clear ii;

end



function F = functional(P,param,f,bound)

alpha = param(1);
beta = param(2);
gamma = param(3);
% F = P/sum(P) - [abs(f).^alpha.*exp(-f.^2/(2*beta^2))]./sum([abs(f).^alpha.*exp(-f.^2/(2*beta^2))]);

P = P/sum(P);

F = [abs(f).^-alpha.*(1-exp(-abs(f).^gamma/(2*beta^gamma)))];
F(~isfinite(F)) = 0;

F = F/sum(F);

F = P-F;

if(~isempty(bound))
    F = F( (f>=bound(1) & f<= bound(2)) | (f>=bound(3) & f<=bound(4))); 
end

F = sum(abs(F));

end

function [c,ceq] = mycon(param,f)

alpha = param(1);
beta = param(2);
gamma = param(3);

c = 0;
P = [abs(f).^-alpha.*(1-exp(-abs(f).^gamma/(2*beta^gamma)))];
P(~isfinite(P)) = 0;

ceq = 1-sum(P);

end

function P = pink_spectrum(param,f)

alpha = param(1);
beta = param(2);
gamma = param(3);
% P = [abs(f).^alpha.*exp(-f.^2/(2*beta^2))]./sum([abs(f).^alpha.*exp(-f.^2/(2*beta^2))]);
P = [abs(f).^-alpha.*(1-exp(-abs(f).^gamma/(2*beta^gamma)))];
P(~isfinite(P)) = 0;

P = P/sum(P);

end
