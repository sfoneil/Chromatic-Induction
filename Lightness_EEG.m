% Here's how the chromatic induction program should be structured:
% The experimental manipulation is to change the background L vs M balance across maybe 16 levels, with 8 more L than the achromatic point and 8 less. We could make the number of levels an input variable.
% On each side of the screen, there should be an achromatic disk. Both disks should be the same size.
% Each of the disks is surrounded by an annulus. One of the annuli is narrower than the other.
% The observer adjusts the L vs M proportion of one disk using the jog wheel on the response device.
% On half of the trials, the left disk is adjusted; on the other half the right disk is adjusted.
% Of the left side adjustment trials, the disk on that side is surrounded by the thinner disk. On the other half, the adjusted disk is surrounded by the thicker disk. And same for the right side adjustment trials.
% We might want to sound a high or low tone to indicate to the participant which side should be adjusted.
% That's basically it, but I don't know how to change the L/M balance while keeping the luminance fixed, which we need to do across both the different background levels and for the disk adjustment. Can talk about that on Friday.
% Let me know if that doesn't make sense or if you have questions.

clear
close all
sca
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);
% white = WhiteIndex(screenNumber);
% grey = white / 2;
% black = BlackIndex(screenNumber);

%addpath('G:\SO Backup\_Post-grad Research\cc');
thisDir = pwd;
%ccDir = 'G:\SO Backup\_Post-grad Research\cc';
ccDir = 'G:\SO Backup\_Post-grad Research\cc'; 
ccDir = 'E:\Experiments\Chromatic Induction';
resultsDir = fullfile(pwd, 'RResults/');
nResults = size(dir(fullfile('resultsDir', '*.xlsx')),1);
%addpath(ccDir)
cd(ccDir);
cc_startup;
cd(thisDir);
load('EMM302_Biosemi_DPP.mat');
%load('SO_home_256_right.mat');

eeg = 0;
debug = 1;
if debug
    screenNum = max(Screen('Screens'));
else
    screenNum = 0;
end

subj = input('Enter subject ID:   ', 's');
d = datetime('now');
subjID = sprintf(d.Year + "-" + d.Month + "-" + d.Day + "_" + ...
    d.Hour + "h" + d.Minute + "m" + ...
    "_" + subj + "_lum");
% cond = input('Enter condition: d=disk, a=annulus', 's');
% if ~any(strcmpi(cond, {'d','a'}))
%     error('Enter d or a!')
% end

%% Parameter defaults
% Monitor
resolution = [1920 1080];
screenSize_cm = [69.7 38.5 79.5]; % l x w x diag
diagLength_cm = sqrt(resolution(1).^2 + resolution(2).^2);
px_per_cm = diagLength_cm / screenSize_cm(3);
viewDistance = 80;
degPerPx(1) = 2*atand((screenSize_cm(1)/resolution(1))/(2*viewDistance));
degPerPx(2) = 2*atand((screenSize_cm(2)/resolution(2))/(2*viewDistance));
degPerPx(3) = 2*atand((screenSize_cm(3)/diagLength_cm)/(2*viewDistance));
pxPerDeg = 1./degPerPx;

% Luminances
lummult = 1;
backgroundLum = 10.^(-1).*lummult; % 51 appears to be max with contrast range of [-80 80]
dotLum = 10.^1 .* lummult;
nLums = 6;
%ringLums = linspace(10^(-0.326), 10^0.201, nLums).*10; % Max appears to be 94 for EMM302 Display++
ringLums = 10.^(linspace(-0.326, 0.402, 6)) .*lummult;
nReps = 6; 
nConds = 2; % 6 lums x 6 reps x 2 conds = 72 trials

backgroundRGB = mw2rgb(trival('MW', [0 0 backgroundLum]));
backgroundRGB = backgroundRGB.Value;
dotRGB = mw2rgb(trival('MW', [0 0 dotLum]));
dotRGB = dotRGB.Value;

% Sizes
dotSize = 2.*[0.35].*pxPerDeg(3); % 1.06 degrees
dotRect = [0 0 dotSize dotSize];
ringSize = 2.*[0.7].*pxPerDeg(3);
ringRect = [0 0 ringSize ringSize];

% Times
stimTime = 2;
isiTime = 1;

% Create trivals
mwDot = trival('MW', [0 0 dotLum; 0 0 dotLum]);
%mwGrayDot = trival('MW', [0 0 dotLum]);
%mwRing = trival('MW', [0 0 ringLum; 0 0 ringLum]);
dotContrastRange = 40; % Randomize range for dot

% Positions
viewSpan = 12./2; % Degrees
ringPos = [resolution(1)./2 - viewSpan.*pxPerDeg(3), ...
    resolution(2)./2; ...
    resolution(1)./2 + viewSpan.*pxPerDeg(3), ...
    resolution(2)./2]';
[xC, yC] = RectCenter([0 0 resolution]);
ringPos = CenterRectOnPoint(ringRect, xC, yC);
dotPos = CenterRectOnPoint(dotRect, xC, yC);
% ringPos = [1920*0.39, 1080*0.5;
%     1920*0.61, 1080*0.5]';

%% Create derived colors
% bkgdTri = [backgroundContrasts', ...
%     zeros(size(backgroundContrasts))', ...
%     ones(size(backgroundContrasts))' .* backgroundLum];
% mwBackground = trival('MW', bkgdTri);
% [rgbBackground, ff] = mw2rgb(mwBackground);
% if any(ff,'all')
%     sca
%     error('One or more background luminances are oversaturated. Lower luminance.')
% end
% rgbGrayDot = mw2rgb(mwGrayDot);
% [rgbRing, ff_ring] = mw2rgb(mwRing);
% if any(ff_ring,'all')
%     sca
%     error('Ring luminances are oversaturated. Lower luminance.')
% end

%% Setup

% Counterbalance order based on number of runs
if mod(nResults, 2) == 0
    mult1 = 1;
    mult2 = 2;
else
    mult1 = 2;
    mult2 = 1;
end

% Create condition blocks
t = table;
for tt = 1:nReps
    % Create first condition
    t1 = table(ringLums', mult1 .* ones(nReps, 1));
    % Shuffle
    t1 = t1(randperm(nReps), :);
    % Create second condition
    t2 = table(ringLums', mult2 .* ones(nReps, 1));
    % Shuffle
    t2 = t2(randperm(nReps), :);
    t = [t; t1; t2];
end
mwLums = trival('MW', [zeros(size(t, 1),1), zeros(size(t, 1),1), ...
    table2array(t(:,1))]);
[rgbLums, ffLums] = mw2rgb(mwLums);
if any(ffLums,'all')
    sca
    error('Ring luminances are oversaturated. Lower luminance.')
end
t(:,3) = num2cell(rgbLums.Value, 2);
t.Properties.VariableNames = {'RingLum', 'Instr', 'RGB'};

escapeKey = KbName('ESCAPE');
% leftKey = KbName('LeftArrow');
% rightKey = KbName('RightArrow');
% setKey = KbName('Space');

%% Instructions
instrDisc = 'THINK about the brightness of each CENTER DISC.\nPress SPACE to continue.';
instrRing = 'THINK about the brightness of the DISC/RING CONTRAST.\nPress SPACE to continue.';


%% Experiment
[win, rect] = PsychImaging('OpenWindow', screenNum, [0.5 0.5 0.5]);
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Flip', win);
ifi = Screen('GetFlipInterval', win);
stimFrames = round(stimTime / ifi);
isiFrames = round(isiTime / ifi);

Screen('TextSize', win, 40);
topPriorityLevel = MaxPriority(win);
if ~debug
    Priority(topPriorityLevel);
    HideCursor;
end

% DrawFormattedText(win, instr, 'center', 'center', [1 1 1]);
% Screen('Flip', win);
% KbWait;
% WaitSecs(0.5);

% Get time stamp
vbl = Screen('Flip', win);

% Experiment block
for tr = 1:size(t,1)
    % Display instructions if 6th trial
    instr = '';
    if mod(tr, nReps)-1 == 0
        if t{tr, 2} == 1
            instr = instrDisc;
        elseif t{tr,2} == 2
            instr = instrRing;
        else
            error('Instruction line error.');
        end
    Screen('FillRect', win, backgroundRGB);
    DrawFormattedText(win, instr, 100, yC, [1 1 1]);
    Screen('Flip', win);
    KbWait;
    WaitSecs(3);
    beep;

    % Quit if ESC
    [keyIsDown,secs, keyCode] = KbCheck(-1);
    if keyCode(escapeKey)
        cd(pwd);
        ShowCursor;
        sca;
        return
    end

    end
    % Get current trial info
    currTrial = t(tr,:);
    currRGB = currTrial{1,3};

        % Draw background
        Screen('FillRect', win, backgroundRGB);

        % Draw rings and dots
        Screen('DrawDots', win, [xC; yC], ringSize, currRGB, [], 1);
        Screen('DrawDots', win, [xC; yC], dotSize, dotRGB, [], 1);

        if debug
            % DrawFormattedText(win, sprintf('Adjust %d\nBkgd = %.4f\nContrast = %.4f', ...
            %     currTrial.Side,currTrial.BackContrast, currContrast), 'center', rect(4)/4, [1 1 1]);
        end
        % Show stimuli
        vbl = Screen('Flip', win);

        % Clear after time is up
        Screen('FillRect', win, backgroundRGB);
        vbl = Screen('Flip', win, vbl + (stimFrames - 0.5) * ifi);

        % ISI
        Screen('FillRect', win, backgroundRGB);
        vbl = Screen('Flip', win, vbl + (isiFrames - 0.5) * ifi);
        
        
end

%% Output and analyze
% Output raw table
%writetable(t, subjID + ".xlsx");

% tAdjSmall = [table2array(t(t.Side==-1 & t.SmallRingIdx==1, ["BackContrast", "Setting"])), ...
%     table2array(t(t.Side==1 & t.SmallRingIdx==2, ["Setting"]))];
% tAdjLarge = [t(t.Side==-1 & t.SmallRingIdx==2, ["BackContrast", "Setting"]), ...
%     t(t.Side==1 & t.SmallRingIdx==1, ["BackContrast", "Setting"])];

% Get small and large adjusts, regardless of side
tAdjSmall = sortrows(t(t.Side==-1 & t.SmallRingIdx==1, ["BackContrast", "Setting"]));
tAdjSmall.Properties.VariableNames = ["BackContrast", "Setting1"];
tAdjSmall = join(tAdjSmall, sortrows(t(t.Side==1 & t.SmallRingIdx==2, ["BackContrast", "Setting"])));

tAdjLarge = sortrows(t(t.Side==-1 & t.SmallRingIdx==2, ["BackContrast", "Setting"]));
tAdjLarge.Properties.VariableNames = ["BackContrast", "Setting1"];
tAdjLarge = join(tAdjLarge, sortrows(t(t.Side==1 & t.SmallRingIdx==1, ["BackContrast", "Setting"])));

% Get descriptives
tAvgSmall = table(tAdjSmall.BackContrast, mean([tAdjSmall.Setting1, tAdjSmall.Setting], 2), ...
    'VariableNames', ["BackContrast", "Small"]);
tAvgLarge = table(tAdjLarge.BackContrast, mean([tAdjLarge.Setting1, tAdjLarge.Setting], 2), ...
    'VariableNames', ["BackContrast", "Large"]);

tSDSmall = table(tAdjSmall.BackContrast, std([tAdjSmall.Setting1, tAdjSmall.Setting], 0, 2), ...
    'VariableNames', ["BackContrast", "Small"]);
tSDLarge = table(tAdjLarge.BackContrast, std([tAdjLarge.Setting1, tAdjLarge.Setting], 0, 2), ...
    'VariableNames', ["BackContrast", "Large"]);

% tAdjSmall = [t(t.Side==-1 & t.SmallRingIdx==1, ["BackContrast", "Setting"]), ...
%     t(t.Side==1 & t.SmallRingIdx==2, "Setting")];
% tAdjLarge = [t(t.Side==-1 & t.SmallRingIdx==2, ["BackContrast", "Setting"]), ...
%     t(t.Side==1 & t.SmallRingIdx==1, "Setting")];

% Average by side
%
% tAdjSmall = sortrows(tAdjSmall);
% tAdjLarge = sortrows(tAdjLarge);
% tAdjSmall.Properties.VariableNames = ["BackContrast", "SmallSet"];
% tAdjLarge.Properties.VariableNames = ["BackContrast", "LargeSet"];

% Compute descriptives
%avg = mean([table2array(tAdjSmall(:,2)), table2array(tAdjLarge(:,2))], 2);
%sd = std([table2array(tAdjSmall(:,2)), table2array(tAdjLarge(:,2))], 0, 2);

% Output to Excel
writetable(t, subjID + ".xlsx")
% writetable(sortrows(t(:, ["BackContrast", "Side", "SmallRingIdx", "Setting"])), subjID + ".xlsx", 'Range', "A1:D" + string(size(t,1)+1));
% writetable(tAvgSmall, subjID + ".xlsx", 'Range', "F1:G" + size(tAvgSmall,1)+1);
% writetable(tSDSmall(:,2), subjID + ".xlsx", 'Range', "H1:H" + size(tSDSmall,1)+1);
% writetable(tAvgLarge(:,2), subjID + ".xlsx", 'Range', "I1:I" + size(tAvgLarge,1)+1);
% writetable(tSDLarge(:,2), subjID + ".xlsx", 'Range', "J1:J" + size(tSDLarge,1)+1);

% tOut = [tAdjSmall, tAdjLarge(:,2)];
% writetable(tOut, subjID + ".xlsx", 'Range', "A1:C" + string(size(tOut,1)+1));
% writecell({'Mean', 'SD'}, subjID + ".xlsx", 'Range', "E1:F1");
% writematrix([avg, sd], subjID + ".xlsx", 'Range', "E2:F"+ string(size(avg,1)+1));

ShowCursor;
sca
