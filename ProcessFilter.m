
function OutputFile = ProcessFilter(sProcess, sInput)
OutputFile = [];
fileTag = '';

% ===== SELECT CHANNELS =====
% Read the channel file
if ~isempty(sInput.ChannelFile)
    ChannelMat = in_bst_channel(sInput.ChannelFile);
else
    ChannelMat = [];
end
% Specific selection
if ismember(sInput.FileType, {'data', 'raw'}) && ~isempty(sInput.ChannelFile) && isfield(sProcess.options, 'sensortypes') && ~isempty(sProcess.options.sensortypes)
    % Get channel indices
    iSelRows = channel_find(ChannelMat.Channel, sProcess.options.sensortypes.Value);
    % If no selection: file not processed
    if isempty(iSelRows)
        bst_report('Error', sProcess, sInput, ['Selected sensor types are not available in file "' sInput.FileName '".']);
        return;
    end
    AllSensorTypes = unique({ChannelMat.Channel(iSelRows).Type});
    % All the signals
else
    iSelRows = [];
    AllSensorTypes = [];
end

% ===== LOAD FILE =====
% Raw file: do not load full file
if strcmpi(sInput.FileType, 'raw')
    isLoadFull = 0;
else
    isLoadFull = 1;
end
% Read input files
[sMat, matName] = in_bst(sInput.FileName, [], isLoadFull);

sInput.Measure = [];
% Is this a continuous file?
isRaw = isstruct(sMat.(matName));

% Absolute values of sources / norm or unconstrained sources

% Get data matrix
matValues = sMat.(matName);
% Get std matrix
if isfield(sMat, 'Std') && ~isempty(sMat.Std)
    stdValues = sMat.Std;
else
    stdValues = [];
end

% Progress bar comment
txtProgress = ['Running process: ' sProcess.Comment '...'];
% Copy channel flag information
if isfield(sMat, 'ChannelFlag')
    sInput.ChannelFlag = sMat.ChannelFlag;
end
% Copy nAvg information
if isfield(sMat, 'nAvg') && ~isempty(sMat.nAvg)
    sInput.nAvg = sMat.nAvg;
else
    sInput.nAvg = 1;
end
% Raw files
isReadAll = isRaw && isfield(sProcess.options, 'read_all') && isfield(sProcess.options.read_all, 'Value') && isequal(sProcess.options.read_all.Value, 1);

sFileIn = matValues;
clear matValues;
iEpoch = 1;
nEpochs = length(sFileIn.epochs);
% Get size of input data
nRow = length(sMat.ChannelFlag);
nCol = length(sMat.Time);
% Get subject
sSubject = bst_get('Subject', sInput.SubjectName);
% ERROR: File does not exist
if ~file_exist(sFileIn.filename)
    bst_report('Error', sProcess, sInput, [...
        'This file has been moved, deleted, is used by another program,' 10 ...
        'or is on a drive that is currently not connected to your computer.']);
    return;
end

% ERROR: SSP cannot be applied for channel/channel processing
if ismember(1,sProcess.processDim) && ~isReadAll && ~isempty(ChannelMat.Projector) && any([ChannelMat.Projector.Status] == 1)
    bst_report('Error', sProcess, sInput, [...
        'This file contains SSP projectors, which require all the channels to be read at the same time.' 10 ...
        'To process this file, you have the following options: ' 10 ...
        '  1) Check the option "Process the entire file at once" (possible only if the entire file fits in memory).' 10 ...
        '  2) Run the process "Artifacts > Apply SSP & CTF compensation" first to save a compensated file.' 10 ...
        '  3) Delete the SSP from this file, process it, then recalculate the SSP on the new file.']);
    return;
end

% Prepare import options
% NOTE: FORCE READING CLEAN DATA (CTF compensators + Previous SSP)
ImportOptions = db_template('ImportOptions');
ImportOptions.ImportMode      = 'Time';
ImportOptions.DisplayMessages = 0;
ImportOptions.UseCtfComp      = 0;
ImportOptions.UseSsp          = 0;
ImportOptions.RemoveBaseline  = 'no';
% Force reading of the entire RAW file at once
if isReadAll
    bst_progress('text', [txtProgress, ' [reading]']);
    FullFileMat = in_fread(sFileIn, ChannelMat, iEpoch, [], [], ImportOptions);
end
% If native file with multiple epochs: ERROR
if isRaw && (nEpochs > 1)
    bst_report('Error', sProcess, sInput, 'Impossible to process native epoched/averaged files. Please import them in database or convert them to continuous.');
    return;
end
%     % Build output file tag
fileTag = [fileTag, '_', strtrim(strrep(sProcess.FileTag,'|',''))];
fileTag = getTag(fileTag);
% Get file type
fileType = file_gettype(sInput.FileName);

% ===== OVERWRITE ? =====
isOverwrite = isfield(sProcess.options, 'overwrite') && sProcess.options.overwrite.Value;
% Overwrite required: check if it is doable
if isOverwrite
    % Ignore overwrite for links
    if strcmpi(fileType, 'link')
        isOverwrite = 0;
        bst_report('Warning', sProcess, sInput, 'Cannot overwrite links.');
    end
end

% ===== OUTPUT FILENAME =====
% Protocol folders and processing options
ProtocolInfo = bst_get('ProtocolInfo');
ProcessOptions = bst_get('ProcessOptions');
% If file is a raw link: create new condition
% Get input raw path and name
[rawPathIn, rawBaseIn] = bst_fileparts(sFileIn.filename);
% Make sure that there are not weird characters in the folder names
rawBaseIn = file_standardize(rawBaseIn);
% Save in original folder

if isfield(sFileIn, 'condition') && ~isempty(sFileIn.condition)
    newCondition = ['@raw', sFileIn.condition, '_preproc'];
else
    newCondition = ['@raw', rawBaseIn, '_preproc'];
end
% Get new condition name
%     newStudyPath = file_unique(bst_fullfile(ProtocolInfo.STUDIES, sInput.SubjectName, newCondition));
newStudyPath = (bst_fullfile(ProtocolInfo.STUDIES, sInput.SubjectName, newCondition));

isLinkRaw = 0;

if(strcmp(sFileIn.format,'BST-BIN'))
    rawBaseOut = rawBaseIn;
else
    % Output file name derives from the condition name
    [tmp, rawBaseOut] = bst_fileparts(newStudyPath);
    isLinkRaw = 1;
end

rawBaseOut = strrep(rawBaseOut, '@raw', '');
% Full output filename
RawFileOut = bst_fullfile(newStudyPath, [rawBaseOut fileTag '.bst']);
if(not(isOverwrite))
    RawFileOut = file_unique(RawFileOut);
end
RawFileFormat = 'BST-BIN';

% Get new condition name
[tmp, ConditionName] = bst_fileparts(newStudyPath, 1);
[sOutputStudy,iOutputStudy] = bst_get('Study',ConditionName);

if(isempty(iOutputStudy))% Create output condition
    iOutputStudy = db_add_condition(sInput.SubjectName, ConditionName);
    if isempty(iOutputStudy)
        bst_report('Error', sProcess, sInput, ['Output folder could not be created:' 10 newPath]);
        return;
    end
    sOutputStudy = bst_get('Study', iOutputStudy);
end
% Get output study
% Full file name
[dum,file_suffix] = fileparts(RawFileOut);
MatFile = bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sOutputStudy.FileName), ['data_0raw_' file_suffix '.mat']);

% ===== SPLIT IN BLOCKS =====
OutMeasure = [];
OutputMat = [];
OutputStd = [];
% Get maximum size of a data block
MaxSize = ProcessOptions.MaxBlockSize;
% Split the block size in rows and columns
if (nRow * nCol > MaxSize) && ~isempty(sProcess.processDim)
    % Split max block by row blocks
    if ismember(1, sProcess.processDim)
        % Split by row and col blocks
        if (nCol > MaxSize) && ismember(2, sProcess.processDim)
            BlockSizeRow = 1;
            BlockSizeCol = MaxSize;
            % Split only by row blocks
        else
            BlockSizeRow = max(floor(MaxSize / nCol), 1);
            BlockSizeCol = nCol;
        end
        % Split max block by col blocks
    elseif ismember(2, sProcess.processDim)
        BlockSizeRow = nRow;
        BlockSizeCol = max(floor(MaxSize / nRow), 1);
    end
    % Adapt block size to FIF block size
    if (BlockSizeCol < nCol) && isRaw && strcmpi(sFileIn.format, 'FIF') && isfield(sFileIn.header, 'raw') && isfield(sFileIn.header.raw, 'rawdir') && ~isempty(sFileIn.header.raw.rawdir)
        fifBlockSize = double(sFileIn.header.raw.rawdir(1).nsamp);
        BlockSizeCol = fifBlockSize * max(1, round(BlockSizeCol / fifBlockSize));
    end
else
    BlockSizeRow = nRow;
    BlockSizeCol = nCol;
end
% Split data in blocks
nBlockRow = ceil(nRow / BlockSizeRow);
nBlockCol = ceil(nCol / BlockSizeCol);
% Get current progress bar position
progressPos = bst_progress('get');
prevPos = 0;
% Display console message
if (nBlockRow > 1) && (nBlockCol > 1)
    disp(sprintf('BST> %s: Processing %d blocks of %d signals and %d time points.', sProcess.Comment, nBlockCol * nBlockRow, BlockSizeRow, BlockSizeCol));
elseif (nBlockRow > 1)
    disp(sprintf('BST> %s: Processing %d blocks of %d signals.', sProcess.Comment, nBlockRow, BlockSizeRow));
elseif (nBlockCol > 1)
    disp(sprintf('BST> %s: Processing %d blocks of %d time points.', sProcess.Comment, nBlockCol, BlockSizeCol));
end

% ===== PROCESS BLOCKS =====
isFirstLoop = 1;
% Loop on row blocks
for iBlockRow = 1:nBlockRow
    % Indices of rows to process
    iRow = 1 + (((iBlockRow-1)*BlockSizeRow) : min(iBlockRow * BlockSizeRow - 1, nRow - 1));
    % Process only the required rows
    if ~isempty(iSelRows)
        [tmp__, iRowProcess] = intersect(iRow, iSelRows);
    end
    % Loop on col blocks
    for iBlockCol = 1:nBlockCol
        % Indices of columns to process
        iCol = 1 + (((iBlockCol-1)*BlockSizeCol) : min(iBlockCol * BlockSizeCol - 1, nCol - 1));
        % Progress bar
        newPos = progressPos + round(((iBlockRow - 1) * nBlockCol + iBlockCol) / (nBlockRow * nBlockCol) * 100);
        if (newPos ~= prevPos)
            bst_progress('set', newPos);
            prevPos = newPos;
        end
        
        % === GET DATA ===
        % Read values
        bst_progress('text', [txtProgress, ' [reading]']);
        % Read block
        if isReadAll
            sInput.A = FullFileMat(iRow, iCol);
        else
            SamplesBounds = sFileIn.prop.samples(1) + iCol([1,end]) - 1;
            sInput.A = in_fread(sFileIn, ChannelMat, iEpoch, SamplesBounds, iRow, ImportOptions);
        end
        sInput.Std = [];
        % Progress bar: processing
        bst_progress('text', [txtProgress, ' [processing]']);
        % Set time vector in input
        sInput.TimeVector = sMat.Time(iCol);
        
        % === PROCESS ===
        % Send indices to the process
        sInput.iBlockRow = iBlockRow;
        sInput.iBlockCol = iBlockCol;
        % Process all rows
        if isempty(iSelRows) || isequal(iRowProcess, 1:size(sInput.A,1))
            sInput.iRowProcess = (1:size(sInput.A,1))';
            sInput = sProcess.Function('NoCatch', 'Run', sProcess, sInput);
            % Process only a subset of rows
        elseif ~isempty(iRowProcess)
            sInput.iRowProcess = iRowProcess;
            tmp1 = sInput.A;
            % Main data matrix
            sInput.A = sInput.A(iRowProcess,:,:);
            % Standard error
            if ~isempty(sInput.Std)
                tmp2 = sInput.Std;
                sInput.Std = sInput.Std(iRowProcess,:,:);
            end
            % Process file
            sInput = sProcess.Function('NoCatch', 'Run', sProcess, sInput);
            % Get results
            if ~isempty(sInput)
                % Main data matrix
                tmp1(iRowProcess,:,:) = sInput.A;
                sInput.A = tmp1;
                % Standard error
                if ~isempty(sInput.Std)
                    tmp2(iRowProcess,:,:) = sInput.Std;
                    sInput.A = tmp2;
                end
            end
        end
        
        % If an error occured
        if isempty(sInput)
            return;
        end
        
        % === INITIALIZE OUTPUT ===
        % Split along columns (time): No support for change in sample numbers (resample)
        if ismember(2, sProcess.processDim)
            nOutTime = nCol;
            iOutTime = iCol;
            % All the other options (split by row, no split): support for resampling
        else
            nOutTime = length(sInput.TimeVector);
            iOutTime = iCol(1) - 1 + (1:length(sInput.TimeVector));
        end
        
        % Create output variable
        if isFirstLoop
            isFirstLoop = 0;
            bst_progress('text', [txtProgress, ' [creating new file]']);
            % Output measure
            if isfield(sInput, 'Measure')
                OutMeasure = sInput.Measure;
            end
            % RAW: Create a new raw file to store the results
            % Create an empty Brainstorm-binary file
            [sFileOut, errMsg] = out_fopen(RawFileOut, RawFileFormat, sFileIn, ChannelMat);
            % Error processing
            if isempty(sFileOut) && ~isempty(errMsg)
                bst_report('Error', sProcess, sInput, errMsg);
                return;
            elseif ~isempty(errMsg)
                bst_report('Warning', sProcess, sInput, errMsg);
            end
            % Did time definition change?
            isTimeChange = ~ismember(2, sProcess.processDim) && ~isequal(sInput.TimeVector, sMat.Time) && (isRaw || ~((size(matValues,2) == 1) && (length(sMat.Time) == 2)));
            % Output time vector
            if isTimeChange
                
                OldFreq = 1./(sMat.Time(2) - sMat.Time(1));
                NewFreq = 1./(sInput.TimeVector(2) - sInput.TimeVector(1));
                if(NewFreq ~= OldFreq)
                    sFileOut.prop.sfreq = NewFreq;
                end
                sFileOut.prop.times = [sInput.TimeVector(1), sInput.TimeVector(end)];
                sFileOut.prop.samples = round(sFileOut.prop.times*sFileOut.prop.sfreq);
%                 sFileOut.header.sfreq = NewFreq;
%                 sFileOut.header.nsamples =  length(sInput.TimeVector);
%                 sFileOut.header.epochsize =  NewFreq;
                % If there are events: update the time and sample indices
                if isfield(sFileIn, 'events') && ~isempty(sFileIn.events)
                    OldFreq = 1./(sMat.Time(2) - sMat.Time(1));
                    if (length(sInput.TimeVector) >= 2)
                        sFileOut.events = panel_record('ChangeTimeVector', sFileIn.events, OldFreq, sInput.TimeVector);
                    else
                        sFileOut.events = [];
                    end
                end
                % Save new time vector
                OutTime = sInput.TimeVector;
            else
                OutTime = sMat.Time;
            end
            % Output channel file
            ChannelMatOut = ChannelMat;
        end
        
        % === SAVE VALUES ===
        bst_progress('text', [txtProgress, ' [writing]']);
        if isReadAll
            FullFileMat(iRow, iCol) = sInput.A;
        else
            % Indices to write
            SamplesBounds = sFileOut.prop.samples(1) + iOutTime([1,end]) - 1;
            % Write block
            sFileOut = out_fwrite(sFileOut, ChannelMatOut, iEpoch, SamplesBounds, iRow, sInput.A);
        end
    end
end

% Save all the RAW file at once
if isReadAll
    sFileOut = out_fwrite(sFileOut, ChannelMatOut, iEpoch, [], [], FullFileMat);
end

% ===== CREATE OUTPUT STRUCTURE =====
% If there is a DataFile link, and the time definition changed, and results is not static: remove link
if isfield(sMat, 'DataFile') && ~isempty(sMat.DataFile)
    if ~isequal(sMat.Time, OutTime) && (length(OutTime) > 2)
        sMat.DataFile = [];
    else
        sMat.DataFile = file_short(sMat.DataFile);
    end
end
% Output time vector
sMat.Time = [OutTime(1), OutTime(end)];
% Output measure
if ~isempty(OutMeasure)
    sMat.Measure = OutMeasure;
end
% Set data fields
% Remove the string: "Link to raw file"
sMat.Comment = strrep(sMat.Comment, 'Link to raw file', 'Raw');
% sMat.Time = [sMat.Time(1), sMat.Time(end)];
sMat.F = sFileOut;

% Comment: forced in the options
if isfield(sProcess.options, 'Comment') && isfield(sProcess.options.Comment, 'Value') && ~isempty(sProcess.options.Comment.Value)
    sMat.Comment = sProcess.options.Comment.Value;
else
    % Add file tag
    if isfield(sInput, 'FileTag') && ~isempty(sInput.FileTag)
        sMat.Comment = [sMat.Comment, ' ', sInput.FileTag];
    else
        sMat.Comment = [sMat.Comment, ' ', sProcess.FileTag];
    end
end
% If data + changed data type
if isfield(sInput, 'DataType') && ~isempty(sInput.DataType) && isfield(sMat, 'DataType')
    sMat.DataType = sInput.DataType;
end
if isfield(sInput, 'ColormapType') && ~isempty(sInput.ColormapType)
    sMat.ColormapType = sInput.ColormapType;
end
if isfield(sInput, 'Function') && ~isempty(sInput.Function)
    sMat.Function = sInput.Function;
end
% ChannelFlag
if isfield(sInput, 'ChannelFlag') && ~isempty(sInput.ChannelFlag)
    sMat.ChannelFlag = sInput.ChannelFlag;
end

% History: Process name + options
if isfield(sInput, 'HistoryComment') && ~isempty(sInput.HistoryComment)
    HistoryComment = [func2str(sProcess.Function) ': ' sInput.HistoryComment];
else
    HistoryComment = [func2str(sProcess.Function) ': ' sProcess.Function('FormatComment', sProcess)];
end
sMat = bst_history('add', sMat, 'process', HistoryComment);

% ===== SAVE FILE =====
% Save new file
bst_save(MatFile, sMat, 'v6');
% If no default channel file: create new channel file
if isRaw && (sSubject.UseDefaultChannel == 0)
    db_set_channel(iOutputStudy, ChannelMatOut, 2, 0);
end

% ===== REGISTER IN DATABASE =====
% Register in database
if isOverwrite
    if(isLinkRaw)
        db_reload_studies(iOutputStudy);
    else
        db_add_data(iOutputStudy, MatFile, sMat, sInput.iItem);        
    end
%     
%     [sStudy, iStudy, iData] = bst_get('DataFile',    MatFile, iOutputStudy);
%     if(isempty(iData))
%         db_add_data(iOutputStudy, MatFile, sMat);
%     elseif
%         db_add_data(iOutputStudy, MatFile, sMat, sInput.iItem);        
%     else
%         db_reload_studies(iOutputStudy);
%     end
else
    db_add_data(iOutputStudy, MatFile, sMat);
end
% Return new file
OutputFile = MatFile;
end

function  fileTag = getTag(fileTag)

if(find(ismember({'_bandpass','_notch'},fileTag)))
    fileTag = '_f';
elseif(find(ismember({'_resample'},fileTag)))
    fileTag = '_r';
elseif(find(ismember({'_channel_artefact'},fileTag)))
    fileTag = '_a';
elseif(find(ismember({'_concat','_extract'},fileTag)))
    fileTag = '_t';
end
end
