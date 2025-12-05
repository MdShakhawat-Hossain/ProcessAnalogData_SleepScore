function Generate_Sleep_Figures_Analog(procDataFileID,saveFigs)
%GENERATE_SLEEP_FIGURES_ANALOG  Summary figure for manual sleep inspection.
%
%   Generate_Sleep_Figures_Analog(procDataFileID,saveFigs)
%
%   LAYOUT
%       ax1: Force + binForce (top)
%       ax2: EMG power (middle)
%       ax3: ECoG spectrogram (bottom, log y-axis, with colorbar)
%
%   INPUTS
%       procDataFileID  – path to *ProcData.mat* file
%       saveFigs        – 'y' or 'n' (default 'n')
%
%   Notes
%       • Uses ProcData.notes.dsFs if available, else falls back to 60 Hz
%       • Works with both KECK-style ProcData.* and older ProcData.data.* formats
%       • Loads SpecData* file (SpecData.ECoG.*) and plots spectrogram
%       • Draws hypnogram strips (Awake/NREM/REM) above each subplot
%         using ProcData.sleep.logicals.Manual.* when available
%         (else falls back to ProcData.sleep.behavState)

if nargin < 2 || isempty(saveFigs)
    saveFigs = 'n';
end

%% ------------------------------------------------------------------------
% Load ProcData and determine timing info
% -------------------------------------------------------------------------
S = load(procDataFileID,'ProcData');
ProcData = S.ProcData;

% Sampling rate
Fs = [];
if isfield(ProcData,'notes') && isfield(ProcData.notes,'dsFs')
    Fs = ProcData.notes.dsFs;
end
if isempty(Fs) || ~isnumeric(Fs) || ~isfinite(Fs)
    Fs = 60;   % fallback
end

% Pick a reference signal to determine length
refLen = [];
if isfield(ProcData,'ECoG_DS_norm')
    refLen = numel(ProcData.ECoG_DS_norm);
elseif isfield(ProcData,'ECoG_DS')
    refLen = numel(ProcData.ECoG_DS);
else
    error('Generate_Sleep_Figures_Analog:CannotDetermineDuration', ...
        'Cannot determine trial duration from available ProcData fields.');
end

trialDur = refLen / Fs;
t        = (0:refLen-1)/Fs;

%% ------------------------------------------------------------------------
% Extract force, binForce, EMG power (robust, shared with scoring figure)
% -------------------------------------------------------------------------
% ===== FORCE =====
if isfield(ProcData,'forceSensor_norm')
    force = ProcData.forceSensor_norm(:);
    forceLabel = 'Force (norm)';
elseif isfield(ProcData,'forceSensor')
    force = ProcData.forceSensor(:);
    forceLabel = 'Force (V)';
else
    force = nan(refLen,1);
    forceLabel = 'Force';
end
force = force(:);
if numel(force) > refLen
    force = force(1:refLen);
elseif numel(force) < refLen
    if ~isempty(force)
        force(end+1:refLen) = force(end);
    else
        force = nan(refLen,1);
    end
end

% ===== BIN FORCE =====
if isfield(ProcData,'binForceSensor')
    binForce = ProcData.binForceSensor(:);
else
    binForce = false(refLen,1);
end
binForce = logical(binForce(:));
if numel(binForce) > refLen
    binForce = binForce(1:refLen);
elseif numel(binForce) < refLen
    tmp = false(refLen,1);
    tmp(1:numel(binForce)) = binForce;
    binForce = tmp;
end

% ===== EMG POWER =====
emgPower = [];
emgLabel = 'EMG power';

% Preferred: ProcData.EMG.emgPower_norm / emgPower
if isfield(ProcData,'EMG') && isstruct(ProcData.EMG)
    if isfield(ProcData.EMG,'emgPower_norm')
        emgPower = ProcData.EMG.emgPower_norm(:);
        emgLabel = 'EMG power (norm)';
    elseif isfield(ProcData.EMG,'emgPower')
        emgPower = ProcData.EMG.emgPower(:);
        emgLabel = 'EMG power';
    end
end

if isempty(emgPower)
    emgPower = nan(refLen,1);
    emgLabel = 'EMG power';
end

emgPower = emgPower(:);
if numel(emgPower) > refLen
    emgPower = emgPower(1:refLen);
elseif numel(emgPower) < refLen
    if ~isempty(emgPower)
        emgPower(end+1:refLen) = emgPower(end);
    else
        emgPower = nan(refLen,1);
    end
end

%% ------------------------------------------------------------------------
% Sleep score (behavioral state) for hypnogram (Awake/NREM/REM)
% -------------------------------------------------------------------------
AwakeStage = [];
NREMStage  = [];
REMStage   = [];
SleepDummy = [];
binDur     = 5;   % seconds per bin (as in your original logic)
sleepOK    = false;

if isfield(ProcData,'sleep')
    % --- Preferred: Manual logicals (awake/nrem/rem) ---
    if isfield(ProcData.sleep,'logicals') && ...
       isfield(ProcData.sleep.logicals,'Manual') && ...
       isfield(ProcData.sleep.logicals.Manual,'awakeLogical') && ...
       isfield(ProcData.sleep.logicals.Manual,'nremLogical')  && ...
       isfield(ProcData.sleep.logicals.Manual,'remLogical')

        AwakeStage = double(ProcData.sleep.logicals.Manual.awakeLogical(:));
        NREMStage  = double(ProcData.sleep.logicals.Manual.nremLogical(:));
        REMStage   = double(ProcData.sleep.logicals.Manual.remLogical(:));

        % zero -> NaN (your original pattern)
        AwakeStage(AwakeStage==0) = NaN;
        NREMStage(NREMStage==0)   = NaN;
        REMStage(REMStage==0)     = NaN;

        TableSize  = numel(AwakeStage);
        SleepDummy = 1:binDur:(TableSize*binDur);   % like your original "1:5:trialDuration_sec"
        sleepOK    = true;

    % --- Fallback: behavState labels -> Awake/NREM/REM ---
    elseif isfield(ProcData.sleep,'behavState')
        sleepLabels = ProcData.sleep.behavState;
        if isstring(sleepLabels)
            sleepLabels = cellstr(sleepLabels);
        end
        numBins = numel(sleepLabels);

        AwakeStage = nan(numBins,1);
        NREMStage  = nan(numBins,1);
        REMStage   = nan(numBins,1);

        for b = 1:numBins
            lab = sleepLabels{b};
            if isempty(lab), continue; end
            switch lab
                case {'Not Sleep','Wake','Wakefulness'}
                    AwakeStage(b) = 1;
                case {'NREM Sleep','NREM','Non-REM'}
                    NREMStage(b)  = 1;
                case {'REM Sleep','REM'}
                    REMStage(b)   = 1;
            end
        end

        SleepDummy = 1:binDur:(numBins*binDur);
        sleepOK    = true;
    end
end

% Ensure SleepDummy is row, others column
SleepDummy = SleepDummy(:).';
AwakeStage = AwakeStage(:);
NREMStage  = NREMStage(:);
REMStage   = REMStage(:);

%% ------------------------------------------------------------------------
% Load spectrogram (same style as scoring figure)
% -------------------------------------------------------------------------
Tspec = []; Fspec = []; Sspec = [];
isNormSpec = false;

[folder, baseName, ~] = fileparts(procDataFileID);
prefix = regexprep(baseName,'_ProcData$','');
specFile = '';
cand = dir(fullfile(folder, [prefix '*Spec*.mat']));
if ~isempty(cand)
    specFile = fullfile(folder, cand(1).name);
end

if ~isempty(specFile) && isfile(specFile)
    try
        L = load(specFile);
        % Prefer SpecData.ECoG.normS if present
        if isfield(L,'SpecData') && isfield(L.SpecData,'ECoG')
            SD = L.SpecData;
            if isfield(SD.ECoG,'normS')
                Sspec = SD.ECoG.normS;
                Fspec = SD.ECoG.F;
                Tspec = SD.ECoG.T;
                isNormSpec = true;
            elseif all(isfield(SD.ECoG,{'S','F','T'}))
                Sspec = SD.ECoG.S;
                Fspec = SD.ECoG.F;
                Tspec = SD.ECoG.T;
            end
        end

        % Fallback: any struct with S/F/T
        if isempty(Sspec)
            fn = fieldnames(L);
            for w = 1:numel(fn)
                val = L.(fn{w});
                if isstruct(val) && all(isfield(val,{'S','F','T'}))
                    Sspec = val.S.*100;
                    Fspec = val.F;
                    Tspec = val.T;
                    break
                end
            end
        end
    catch
        Tspec = []; Fspec = []; Sspec = [];
    end
end

%% ------------------------------------------------------------------------
% Create figure and subplots
% -------------------------------------------------------------------------
figHandle = figure('Name','Sleep figure','Color','w');
set(figHandle,'Units','normalized','Position',[0.1 0.1 0.8 0.75]);

ax1 = subplot(4,1,1,'Parent',figHandle);  % Force + binForce
ax2 = subplot(4,1,2,'Parent',figHandle);  % EMG power
ax3 = subplot(4,1,3:4,'Parent',figHandle);  % Spectrogram

%% ------------------------------------------------------------------------
% Top: Force + binForce (keep yyaxis style)
% -------------------------------------------------------------------------
axes(ax1); hold(ax1,'on');

if ~isempty(force)
    yyaxis(ax1,'left');
    plot(ax1,t, force, 'Color',[0.2 0.7 0.2], 'LineWidth',1);
    ylabel(ax1,forceLabel);
end

if ~isempty(binForce)
    idx = find(binForce);
    if ~isempty(idx)
        yyaxis(ax1,'right');
        scatter(ax1, t(idx), ones(size(idx)), 15, 'r', 'filled');
        ylabel(ax1,'binForce');
        yyaxis(ax1,'left');
    end
end

xlim(ax1,[0 trialDur]);
set(ax1,'XTickLabel',[]);
box(ax1,'on');
title(ax1,strrep(procDataFileID,'_','\_'));

%% ------------------------------------------------------------------------
% Middle: EMG power
% -------------------------------------------------------------------------
axes(ax2); hold(ax2,'on');
if ~isempty(emgPower)
    plot(ax2,t, emgPower, 'Color',[0 0 0], 'LineWidth',1);
    ylabel(ax2,emgLabel);
else
    text(ax2,0.5,0.5,'EMG power not available','Units','normalized',...
        'HorizontalAlignment','center','VerticalAlignment','middle',....
        'Color',[0.5 0.5 0.5],'FontAngle','italic');
    ylabel(ax2,'EMG');
end
xlim(ax2,[0 trialDur]);
set(ax2,'XTickLabel',[]);
box(ax2,'on');

%% ------------------------------------------------------------------------
% Bottom: Spectrogram (percent change if normalized)
% -------------------------------------------------------------------------
axes(ax3); hold(ax3,'on');

if ~isempty(Sspec) && ~isempty(Tspec) && ~isempty(Fspec)

    % -----------------------------------------
    % CASE 1: Kevin-style normalized: (S - mu)/mu
    % -----------------------------------------
    if isfield(ProcData,'RestingBaseline') && ...
       isfield(ProcData.RestingBaseline,'spg') && ...
       isfield(ProcData.RestingBaseline.spg,'mean') && ...
       isNormSpec

        % ΔS/S0  (fractional change)
        fracChange = Sspec;                          % F x T

        % Convert to percent
        Splot = 100 * fracChange;                    % percent change

        imagesc(ax3, Tspec, Fspec, Splot);
        axis(ax3,'xy');

        % ----- CHOOSE GOOD CLim for percent change -----
        % Suppress everything below -50%, highlight up to +150%
        lo = -100;           % collapse very low power into dark
        hi = 100;           % typical upper range for spindles/gamma
        set(ax3,'CLim',[lo hi]);

        c = colorbar(ax3,'Location','eastoutside');
        c.Label.String = '\Delta power (% of baseline)';

    % -----------------------------------------
    % CASE 2: No normalization → raw S → convert to dB
    % -----------------------------------------
    else
        Sdb = 10*log10(Sspec + eps);

        imagesc(ax3, Tspec, Fspec, Sdb);
        axis(ax3,'xy');

        % Automatic robust contrast for raw power
        p = prctile(Sdb(:),[5 95]);
        if all(isfinite(p))
            set(ax3,'CLim',p);
        end

        c = colorbar(ax3,'Location','eastoutside');
        c.Label.String = 'Power (dB)';
    end

    % -----------------------------------------
    % Common formatting
    % -----------------------------------------
    set(ax3,'YScale','log');
    yticks(ax3,[4 8 15 30]);
    yticklabels(ax3,{'4','8','15','30'});

    xlabel(ax3,'Time (s)');
    ylabel(ax3,'Freq (Hz)');
    xlim(ax3,[0 trialDur]);
    box(ax3,'on');

else
    % -----------------------------------------
    % No spectrogram case
    % -----------------------------------------
    text(ax3,0.5,0.5,'Spectrogram not available','Units','normalized', ...
         'HorizontalAlignment','center','VerticalAlignment','middle', ...
         'Color',[0.5 0.5 0.5],'FontAngle','italic');
    xlabel(ax3,'Time (s)');
    ylabel(ax3,'Freq (Hz)');
end



%% ------------------------------------------------------------------------
% Link x-axes and make widths identical
% -------------------------------------------------------------------------
linkaxes([ax1,ax2,ax3],'x');
xlim(ax1,[0 trialDur]);

pos1 = get(ax1,'Position');
pos2 = get(ax2,'Position');
pos3 = get(ax3,'Position');

pos2(1) = pos1(1); pos2(3) = pos1(3);
pos3(1) = pos1(1); pos3(3) = pos1(3);

set(ax2,'Position',pos2);
set(ax3,'Position',pos3);

%% ------------------------------------------------------------------------
% X-ticks every 10 s, labels every 50 s on spectrogram axis (match scoring)
% -------------------------------------------------------------------------
tickStep = 10;
maxTick  = floor(trialDur/tickStep)*tickStep;
xt = 0:tickStep:maxTick;
xtLbl = cell(size(xt));
for w = 1:numel(xt)
    if mod(xt(w),50) == 0
        xtLbl{w} = num2str(xt(w));
    else
        xtLbl{w} = '';
    end
end
set(ax3,'XTick',xt,'XTickLabel',xtLbl);

%% ------------------------------------------------------------------------
% Hypnogram strips (Awake/NREM/REM) ABOVE each subplot (no overlay on data)
% -------------------------------------------------------------------------
if sleepOK && ~isempty(SleepDummy)
    stripFrac = 0.1;   % thickness of strip as fraction of y-range

    drawHypnogramStrip(ax1, SleepDummy, AwakeStage, NREMStage, REMStage, binDur, stripFrac);
    drawHypnogramStrip(ax2, SleepDummy, AwakeStage, NREMStage, REMStage, binDur, stripFrac);
    drawHypnogramStrip(ax3, SleepDummy, AwakeStage, NREMStage, REMStage, binDur, stripFrac);
end

%% ------------------------------------------------------------------------
% Optional: save figure
% -------------------------------------------------------------------------
if strcmpi(saveFigs,'y')
    [fp,fn,~] = fileparts(procDataFileID);
    outBase = fullfile(fp, [fn '_SleepFigure']);
    set(figHandle,'PaperPositionMode','auto');
    % print(figHandle,[outBase '.pdf'],'-dpdf','-painters');
    savefig(figHandle,[outBase '.fig']);
    saveas(figHandle,[outBase '.png']);
end

close(figHandle)

%% ========================================================================
% Nested helper: draw hypnogram strip above data in a given axis
% ========================================================================
    function drawHypnogramStrip(ax, SleepDummy, AwakeStage, NREMStage, REMStage, binDur, stripFrac)
        axes(ax); hold(ax,'on');

        % Get current limits and extend Y upward
        yl = ylim(ax);
        yRange = yl(2) - yl(1);
        if yRange <= 0
            yRange = 1;
        end

        delta   = stripFrac * yRange;
        yBottom = yl(2);
        yTop    = yl(2) + delta;
        ylim(ax,[yl(1) yTop]);

        % Build stage vector: 1=Awake, 2=NREM, 3=REM
        nBins = numel(SleepDummy);
        stageVec = nan(nBins,1);
        if ~isempty(AwakeStage)
            stageVec(~isnan(AwakeStage)) = 1;
        end
        if ~isempty(NREMStage)
            stageVec(~isnan(NREMStage))  = 2;
        end
        if ~isempty(REMStage)
            stageVec(~isnan(REMStage))   = 3;
        end

        % Draw colored patches per bin (in seconds)
        for w = 1:nBins
            if isnan(stageVec(w)), continue; end

            switch stageVec(w)
                case 1, c = [0 0 0];                  % Awake (black)
                case 2, c = [0.3010 0.7450 0.9330];   % NREM (blue-ish)
                case 3, c = [1 0 0];                  % REM (red)
                otherwise, continue;
            end

            x0 = SleepDummy(w);
            x1 = x0 + binDur;

            patch('Parent',ax, ...
                  'XData',[x0 x1 x1 x0], ...
                  'YData',[yBottom yBottom yTop yTop], ...
                  'FaceColor',c, ...
                  'EdgeColor','none');
        end
    end

end
