function BaselineNormalization_Analog(varargin)
%BASELINENORMALIZATION_ANALOG
%
% Full pipeline for:
%   • Manual baseline selection (per file, stored in *_BaselineInfo.mat)
%   • Automatic resting-baseline detection (2–20 s, binForce = 0)
%   • Compute ONE baseline mean/std per animal per day
%   • Normalize ProcData signals (force, EMG, ECoG_DS, EMG.emgSignal)
%   • Normalize spectrogram (fractional Δ from resting baseline)
%   • Save normalized signals BACK into existing files (no new files created)
%
% IMPORTANT:
%   → No "_Analog" suffix added to saved variables.
%   → Only the FUNCTION name contains "_Analog".
%

%% Parse options
p = inputParser;
addParameter(p,'ForceReEntry',false,@islogical);  % force re-entry of manual baseline windows
addParameter(p,'Folder',pwd,@ischar);
parse(p,varargin{:});
opts = p.Results;

cd(opts.Folder);
files = dir('*_ProcData.mat');
assert(~isempty(files),'No *_ProcData.mat files found.');

fprintf('\n=== Baseline Normalization (Analog) for %d files ===\n',numel(files));

%% Parse animal + day for each file
nFiles = numel(files);
allNames   = cell(nFiles,1);
animalID   = cell(nFiles,1);
dayCode    = cell(nFiles,1);
groupKey   = cell(nFiles,1);

for i = 1:nFiles
    fn = files(i).name;
    allNames{i} = fn;
    [animalID{i}, dayCode{i}] = parseAnimalDay(fn);
    groupKey{i} = [animalID{i} '_' dayCode{i}];
end

[uniqueKeys,~,keyIdx] = unique(groupKey);

%% LOOP OVER GROUPS (animal + day)
for g = 1:numel(uniqueKeys)
    thisKey   = uniqueKeys{g};
    thisMask  = (keyIdx == g);
    theseFiles = allNames(thisMask);
    thisAnimal = animalID{find(thisMask,1,'first')};
    thisDay    = dayCode{find(thisMask,1,'first')};

    fprintf('\n--------------------------------------------------\n');
    fprintf('Animal %s, Day %s: %d files\n', thisAnimal, thisDay, numel(theseFiles));
    fprintf('Collecting rest samples for baseline...\n');

    % containers for this animal-day
    allRest_force      = [];  % column vector
    allRest_emg        = [];  % EMG power
    allRest_emgSignal  = [];  % raw EMG signal
    allRest_ecog       = [];  % ECoG_DS (may be empty)
    allRest_S          = [];  % spectrogram: F x T_all

    %% PASS 1: collect rest samples across ALL files of this day
    for iF = 1:numel(theseFiles)
        fn = theseFiles{iF};
        baseName = erase(fn,'_ProcData.mat');

        fprintf('Collecting rest samples from %s\n', fn);

        % ---------------------- LOAD PROCDATA ----------------------
        S = load(fn,'ProcData');
        ProcData = S.ProcData;

        if isfield(ProcData.notes,'dsFs')
            Fs = ProcData.notes.dsFs;
        else
            Fs = 60;
        end

        % require forceSensor & EMG & binForceSensor
        if ~isfield(ProcData,'forceSensor') || ...
           ~isfield(ProcData,'EMG') || ...
           ~isfield(ProcData.EMG,'emgPower') || ...
           ~isfield(ProcData,'binForceSensor')
            warning('  Missing required fields in %s. Skipping this file for baseline.', fn);
            continue;
        end

        N = numel(ProcData.forceSensor);
        t = (0:N-1)/Fs;

        % ------------------ MANUAL BASELINE WINDOW ------------------
        baselineInfoFile = [baseName '_BaselineInfo.mat'];
        manualBaseline = [];

        if ~exist(baselineInfoFile,'file') || opts.ForceReEntry
            fprintf(' → Selecting manual baseline window...\n');
            manualBaseline = selectManualBaseline(ProcData, fn, t);
            save(baselineInfoFile,'manualBaseline');
        else
            tmp = load(baselineInfoFile);
            manualBaseline = tmp.manualBaseline;
            fprintf(' → Using stored baseline [%g %g] sec\n',manualBaseline(1),manualBaseline(2));
        end

        % -------------------- AUTOMATIC REST EPOCHS --------------------
        fprintf(' → Detecting rest epochs...\n');
        restIdx = detectRestEpochs(ProcData.binForceSensor(:), Fs);

        mask = t >= manualBaseline(1) & t <= manualBaseline(2);
        idxFinal = restIdx & mask(:);

        if ~any(idxFinal)
            warning('  No rest samples in baseline window for %s. Skipping this file for baseline.', fn);
            continue;
        end

        % ----------------- GATHER REST SAMPLES (COLUMN) -----------------
        fprintf(' → Collecting Force/EMG/ECoG baseline samples...\n');

        % Force
        chunkForce = ProcData.forceSensor(idxFinal);
        allRest_force = [allRest_force; chunkForce(:)]; %#ok<AGROW>

        % EMG power
        chunkEmg = ProcData.EMG.emgPower(idxFinal);
        allRest_emg = [allRest_emg; chunkEmg(:)];       %#ok<AGROW>

        % Raw EMG signal (if present)
        if isfield(ProcData.EMG,'emgSignal')
            chunkEmgSig = ProcData.EMG.emgSignal(idxFinal);
            allRest_emgSignal = [allRest_emgSignal; chunkEmgSig(:)]; %#ok<AGROW>
        end

        % ECoG_DS
        if isfield(ProcData,'ECoG_DS')
            chunkEcog = ProcData.ECoG_DS(idxFinal);
            allRest_ecog = [allRest_ecog; chunkEcog(:)]; %#ok<AGROW>
        end

        % -------------------- LOAD SPECTROGRAM --------------------
        fprintf(' → Collecting spectrogram baseline samples...\n');
        specFile = findSpectrogramFile(fn);
        if isempty(specFile)
            warning('  No spectrogram file found for %s. Skipping spectrogram for this file.', fn);
        else
            S2 = load(specFile,'SpecData');
            SpecData = S2.SpecData;

            if ~isfield(SpecData,'ECoG') || ~isfield(SpecData.ECoG,'S') || ...
               ~isfield(SpecData.ECoG,'T') || ~isfield(SpecData.ECoG,'F')
                warning('  SpecData.ECoG missing S/T/F in %s. Skipping spectrogram for this file.', specFile);
            else
                Sraw  = SpecData.ECoG.S;   % F x T
                Tspec = SpecData.ECoG.T;
                restTimes = t(idxFinal);

                % Map resting samples → closest spectrogram time bins
                if ~isempty(restTimes)
                    specIdx = zeros(size(restTimes));
                    for k = 1:numel(restTimes)
                        [~, ix] = min(abs(Tspec - restTimes(k)));
                        specIdx(k) = ix;
                    end
                    specIdx = unique(specIdx);  % unique columns
                    specIdx(specIdx < 1 | specIdx > size(Sraw,2)) = [];

                    if ~isempty(specIdx)
                        specChunk = Sraw(:, specIdx);  % F x M
                        if isempty(allRest_S)
                            allRest_S = specChunk;
                        else
                            allRest_S = [allRest_S, specChunk]; %#ok<AGROW>
                        end
                    else
                        warning('  No valid spectrogram indices matched rest times in %s.', specFile);
                    end
                end
            end
        end
    end % files in this group

    %% CHECK: any data collected?
    if isempty(allRest_force) || isempty(allRest_emg)
        warning('*** No valid rest samples collected for %s. Skipping this animal-day. ***', thisKey);
        continue;
    end

    fprintf(' → Computing baseline stats for Animal %s, Day %s\n', thisAnimal, thisDay);

    % ----------------- COMPUTE BASELINE MEAN / STD -----------------
    baseline = struct();

    baseline.force.mean = mean(allRest_force,'omitnan');
    baseline.force.std  = std(allRest_force,'omitnan');

    baseline.emg.mean   = mean(allRest_emg,'omitnan');
    baseline.emg.std    = std(allRest_emg,'omitnan');

    % EMG raw signal baseline (if collected)
    if ~isempty(allRest_emgSignal)
        baseline.emgSignal.mean = mean(allRest_emgSignal,'omitnan');
        baseline.emgSignal.std  = std(allRest_emgSignal,'omitnan');
    else
        baseline.emgSignal.mean = NaN;
        baseline.emgSignal.std  = NaN;
    end

    if ~isempty(allRest_ecog)
        baseline.ecog.mean = mean(allRest_ecog,'omitnan');
        baseline.ecog.std  = std(allRest_ecog,'omitnan');
    else
        baseline.ecog.mean = NaN;
        baseline.ecog.std  = NaN;
    end

    if ~isempty(allRest_S)
        % linear-space mean/std (kept for compatibility)
        baseline.spg.mean    = mean(allRest_S, 2, 'omitnan');   % F x 1
        baseline.spg.std     = std(allRest_S, 0, 2, 'omitnan'); % F x 1

        % dB-space baseline (optional, if you want)
        allRest_Sdb          = 10*log10(allRest_S + eps);       % F x T_all
        baseline.spg.mean_db = mean(allRest_Sdb, 2, 'omitnan'); % F x 1
        baseline.spg.std_db  = std(allRest_Sdb, 0, 2, 'omitnan');
    else
        baseline.spg.mean    = [];
        baseline.spg.std     = [];
        baseline.spg.mean_db = [];
        baseline.spg.std_db  = [];
    end

    % Save one baseline file per animal-day
    baselineFile = sprintf('%s_%s_RestingBaseline.mat', thisAnimal, thisDay);
    save(baselineFile,'baseline');
    fprintf(' → Saved baseline for %s to %s\n', thisKey, baselineFile);

    %% PASS 2: apply baseline to each file in this animal-day
    fprintf('Applying baseline normalization to files for %s...\n', thisKey);

    for iF = 1:numel(theseFiles)
        fn = theseFiles{iF};
        baseName = erase(fn,'_ProcData.mat');

        fprintf('  Normalizing %s\n', fn);

        % Load ProcData
        S = load(fn,'ProcData');
        ProcData = S.ProcData;

        if isfield(ProcData.notes,'dsFs')
            Fs = ProcData.notes.dsFs;
        else
            Fs = 60;
        end

        % require force / EMG power
        if ~isfield(ProcData,'forceSensor') || ...
           ~isfield(ProcData,'EMG') || ...
           ~isfield(ProcData.EMG,'emgPower')
            warning('    Missing required fields in %s. Skipping normalization for this file.', fn);
            continue;
        end

        % --------- NORMALIZE PROCDATA CHANNELS ---------
        fprintf('    → Normalizing ProcData...\n');

        % Force
        if baseline.force.std > 0
            ProcData.forceSensor_norm = (ProcData.forceSensor - baseline.force.mean) ./ baseline.force.std;
        else
            ProcData.forceSensor_norm = ProcData.forceSensor * NaN;
        end

        % EMG power
        if baseline.emg.std > 0
            ProcData.EMG.emgPower_norm = (ProcData.EMG.emgPower - baseline.emg.mean) ./ baseline.emg.std;
        else
            ProcData.EMG.emgPower_norm = ProcData.EMG.emgPower * NaN;
        end

        % Raw EMG signal
        if isfield(ProcData.EMG,'emgSignal')
            if isfield(baseline,'emgSignal') && baseline.emgSignal.std > 0
                ProcData.EMG.emgSignal_norm = ...
                    (ProcData.EMG.emgSignal - baseline.emgSignal.mean) ./ baseline.emgSignal.std;
            else
                ProcData.EMG.emgSignal_norm = ProcData.EMG.emgSignal * NaN;
            end
        end

        % ECoG_DS
        if isfield(ProcData,'ECoG_DS') && ~isempty(baseline.ecog.mean) && baseline.ecog.std > 0
            ProcData.ECoG_DS_norm = (ProcData.ECoG_DS - baseline.ecog.mean) ./ baseline.ecog.std;
        elseif isfield(ProcData,'ECoG_DS')
            ProcData.ECoG_DS_norm = ProcData.ECoG_DS * NaN;
        end

        % Attach baseline metadata (animal-day level)
        ProcData.RestingBaseline = baseline;
        ProcData.RestingBaseline.animalID = thisAnimal;
        ProcData.RestingBaseline.dayCode  = thisDay;

        save(fn,'ProcData','-append');
        fprintf('    → Saved normalized ProcData to %s\n', fn);

        % --------- NORMALIZE SPECTROGRAM (fractional Δ from rest) ---------
        specFile = findSpectrogramFile(fn);
        if isempty(specFile)
            fprintf('    → No spectrogram file for %s. Skipping spectrogram normalization.\n', fn);
            continue;
        end

        fprintf('    → Normalizing spectrogram (%s)...\n', specFile);

        S2 = load(specFile,'SpecData');
        SpecData = S2.SpecData;

        if ~isfield(SpecData,'ECoG') || ~isfield(SpecData.ECoG,'S') || ...
           ~isfield(SpecData.ECoG,'F') || ~isfield(SpecData.ECoG,'T')
            warning('    SpecData.ECoG missing S/F/T in %s. Skipping spectrogram.', specFile);
        else
            Sraw = SpecData.ECoG.S;                 % F x T, raw power

            if ~isempty(baseline.spg.mean)
                % fractional change from baseline per frequency
                mu = baseline.spg.mean;             % F x 1

                % Guard against zero / negative baselines
                mu(mu <= 0 | isnan(mu)) = NaN;

                % Broadcast across time
                muMat = repmat(mu, 1, size(Sraw,2));   % F x T

                % ΔS / S0
                SpecData.ECoG.normS = (Sraw - muMat)./muMat;

                % OPTIONAL: also store dB ratio
                SpecData.ECoG.normS_dB = 10*log10(Sraw./muMat);
            else
                SpecData.ECoG.normS = Sraw*NaN;
            end

            save(specFile,'SpecData','-append');
            fprintf('    → Saved normalized spectrogram to %s\n', specFile);
        end

    end % files in group

end % groups

fprintf('\n=== Baseline normalization complete ===\n');
end


%% ==================================================================
%% Helpers inside same file
%% ==================================================================

function manualBaseline = selectManualBaseline(ProcData,Fn,t)
% GUI select baseline start/end (integer seconds)
fs = figure('Name',['Manual baseline: ' Fn],'Color','w');

subplot(3,1,1);
plot(t,ProcData.forceSensor); ylabel('Force'); title('Select baseline start/end'); grid on;

subplot(3,1,2);
if isfield(ProcData,'EMG') && isfield(ProcData.EMG,'emgPower')
    plot(t,ProcData.EMG.emgPower); ylabel('EMG power'); grid on;
else
    plot(t,zeros(size(t))); ylabel('EMG power (none)'); grid on;
end

subplot(3,1,3);
if isfield(ProcData,'ECoG_DS')
    plot(t,ProcData.ECoG_DS); ylabel('ECoG DS');
else
    plot(t,zeros(size(t))); ylabel('ECoG DS (none)');
end
xlabel('Time (s)'); grid on;

disp('Click START then END (integer seconds enforced)');
[x,~] = ginput(2);
manualBaseline = sort(round(x));      % enforce integer seconds
manualBaseline(1) = max(manualBaseline(1), 0);
manualBaseline(2) = max(manualBaseline(2), manualBaseline(1)+1);

close(fs);
end


function restIdx = detectRestEpochs(binForce,Fs)
% Rest = binForce == 0 for 2–20 seconds
isRest = ~logical(binForce(:));

% Find continuous regions
d = diff([0; isRest; 0]);
startIdx = find(d==1);
endIdx   = find(d==-1)-1;

restIdx = false(size(isRest));

for k = 1:numel(startIdx)
    L   = endIdx(k) - startIdx(k) + 1;
    dur = L / Fs;
    if dur >= 2 && dur <= 20
        restIdx(startIdx(k):endIdx(k)) = true;
    end
end
end


function specFile = findSpectrogramFile(procFile)
% Find matching spectrogram file for a given *_ProcData.mat
% Example ProcData:  SGNM000_251203_13_16_1103_ProcData.mat
% SpecData pattern:  SGNM000_251203_13_16_1103_SpecData.mat (or *Spec*.mat)

[folder, base, ~] = fileparts(procFile);
prefix = erase(base,'_ProcData');

cand = dir(fullfile(folder, [prefix '*Spec*.mat']));
if isempty(cand)
    specFile = '';
else
    specFile = fullfile(folder, cand(1).name);
end
end


function [animalID, dayCode] = parseAnimalDay(procFile)
% Parse animal ID and day code from filename
% Expected: AnimalID_YYMMDD_..._ProcData.mat
[~, base, ~] = fileparts(procFile);
parts = strsplit(base,'_');
if numel(parts) < 3
    error('Filename %s does not match expected pattern AnimalID_YYMMDD_...', procFile);
end
animalID = parts{1};
dayCode  = parts{2};  % e.g., '251203'
end
