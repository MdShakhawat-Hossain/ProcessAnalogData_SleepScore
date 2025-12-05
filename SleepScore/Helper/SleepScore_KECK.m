function SleepScore_KECK(opts)
%SLEEPSCORE_KECK  Manual sleep scoring pipeline (analog / Keck-style ProcData).
%
%   SleepScore_KECK() will:
%       1) Find all *_ProcData.mat files in the current folder (TrialFolder)
%       2) Derive NBins from the first file using TrialDurationSeconds / BinSizeSeconds
%       3) Launch the manual scoring GUI to create TrainingData + update ProcData
%
%   SleepScore_KECK(OPT) allows overriding defaults:
%       OPT.BinSizeSeconds  - scalar, seconds per bin (default: 5)
%       OPT.TrialFolder     - folder with *_ProcData.mat (default: pwd)
%
%   OUTPUT FILES (per ProcData file):
%       - *_TrainingData.mat        (contains trainingTable.behavState)
%       - *_ProcData.mat            (updated: ProcData.sleep.behavState, etc. –
%                                    done inside the downstream functions)
%
%   NOTE
%       This is a thin wrapper around CreateTrainingDataSet_ShowScore_Analog.
%
%   © Md Shakhawat Hossain / 2025

clc;

%% ───── OPTIONS ─────
if nargin < 1 || isempty(opts), opts = struct; end
opts = fillDefaults(opts, struct( ...
    'BinSizeSeconds', 5, ...
    'TrialFolder',    pwd));

%% ───── ENV / PATH ─────
startDir = pwd;
c = onCleanup(@() safeCd(startDir)); %#ok<NASGU>

cd(opts.TrialFolder);
fprintf('\nSleepScoreMainScript_Analog: starting in %s\n', pwd);

%% ───── FIND ProcData FILES ─────
procList = dir('*_ProcData.mat');
assert(~isempty(procList), 'No *_ProcData.mat files found in %s', opts.TrialFolder);
procDataFileIDs = char({procList.name}');

% Derive NBins from the FIRST file
firstFile = strtrim(procDataFileIDs(1,:));
tmp = load(firstFile, 'ProcData');
assert(isfield(tmp,'ProcData'), 'First ProcData file does not contain ProcData struct.');

trialDur = getTrialDurationSeconds(tmp.ProcData);  % <-- robust helper
BinSize  = opts.BinSizeSeconds;

NBins = trialDur / BinSize;
if abs(NBins - round(NBins)) > 1e-6
    warning('Trial duration (%.3f s) is not an exact multiple of BinSize (%.3f s). Rounding NBins.', ...
        trialDur, BinSize);
end
NBins = max(1, round(NBins));

fprintf('Found %d ProcData files. BinSize=%.3f s → NBins=%d (trialDur≈%.3f s)\n', ...
    size(procDataFileIDs,1), BinSize, NBins, trialDur);

%% ───── RUN MANUAL SCORING ─────
% CreateTrainingDataSet_SleepAnalog(procDataFileIDs, NBins); % default option
CreateTrainingDataSet_SleepAnalog(procDataFileIDs, NBins, true); % if you want to rescore everything

fprintf('\nSleepScoreMainScript_Analog: completed successfully.\n');

%% --save the sleep score figure ---%
saveFigs = 'y';

for a = 1:size(procDataFileIDs,1)
    procDataFileID      = strtrim(procDataFileIDs(a,:));
    Generate_Sleep_Figures_Analog(procDataFileID,saveFigs);
end
    fprintf('\nSleep Score Single Trial Figure Saved.\n');


end

%% ===================================================================
%%                               HELPERS
%% ===================================================================

function opts = fillDefaults(opts, defaults)
f = fieldnames(defaults);
for k = 1:numel(f)
    if ~isfield(opts, f{k}) || isempty(opts.(f{k}))
        opts.(f{k}) = defaults.(f{k});
    end
end
end

function safeCd(target)
try
    cd(target);
catch
    % silently ignore (e.g., if session closed)
end
end

function trialDur = getTrialDurationSeconds(ProcData)
%GETTRIALDURATIONSECONDS  Robustly infer trial duration (seconds) from ProcData.
%
% Tries, in order:
%   1) ProcData.notes.trialDuration_sec
%   2) ProcData.notes.trialDur_s
%   3) ProcData.notes.duration_sec
%   4) Range of ProcData.data.timeVector
%   5) Length / Fs of a reasonable data channel (force, EMG, etc.)
%
% Errors if nothing works.

trialDur = [];

% -------------------- 1) Notes fields --------------------
if isfield(ProcData,'notes')
    n = ProcData.notes;
    candFields = {'trialDuration_sec','trialDur_s','duration_sec','TrialDur_s','TrialDuration_sec'};
    for i = 1:numel(candFields)
        fn = candFields{i};
        if isfield(n, fn) && ~isempty(n.(fn)) && isnumeric(n.(fn))
            trialDur = double(n.(fn));
            if trialDur > 0
                return;
            end
        end
    end
end

% -------------------- 2) timeVector in data --------------------
if isfield(ProcData,'data') && isfield(ProcData.data,'timeVector')
    tv = ProcData.data.timeVector;
    tv = tv(:);
    if numel(tv) > 1 && all(isfinite(tv))
        trialDur = tv(end) - tv(1);
        if trialDur > 0
            return;
        end
    end
end

% -------------------- 3) Length / Fs from a channel --------------------
Fs = [];
if isfield(ProcData,'notes')
    n = ProcData.notes;
    if isfield(n,'dsFs') && ~isempty(n.dsFs)
        Fs = double(n.dsFs);
    elseif isfield(n,'Fs') && ~isempty(n.Fs)
        Fs = double(n.Fs);
    end
end

L = [];
if isfield(ProcData,'data')
    candChan = { ...
        'forceSensor','Force','ForceSensor', ...
        'EMGPower','EMG_power','EMGpower','EMG', ...
        'GCaMP','GCaMP7','GCaMP_F','analogSignal','signal'};
    for i = 1:numel(candChan)
        nm = candChan{i};
        if isfield(ProcData.data, nm) && ~isempty(ProcData.data.(nm))
            x = ProcData.data.(nm);
            if isnumeric(x) || islogical(x)
                L = numel(x);
                break;
            end
        end
    end
end

if ~isempty(Fs) && ~isempty(L) && Fs > 0
    trialDur = L / Fs;
    if trialDur > 0
        return;
    end
end

% -------------------- If we get here, we failed --------------------
error(['getTrialDurationSeconds: Could not infer trial duration from ProcData. ' ...
       'No usable notes fields (trialDuration_sec, etc.), timeVector, or (length / Fs) found.']);
end
