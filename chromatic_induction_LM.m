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
ccDir = 'E:\Experiments\Chromatic Induction';
%addpath('F:\_Post-grad Research\cc')
cd(ccDir);
cc_startup;
cd(thisDir);
load('EMM302_Biosemi_DPP.mat');
%load('SO_home_256_right.mat');

debug = 0;
if debug
    screenNum = max(Screen('Screens'));
else
    screenNum = 0;
end

subj = input('Enter subject ID:   ', 's');
d = datetime('now');
subjID = sprintf(d.Year + "-" + d.Month + "-" + d.Day + "_" + ...
    d.Hour + "h" + d.Minute + "m" + ...
    "_" + subj + "_LM");

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

% Luminances. 
backgroundLum = 45; % 51 appears to be max with contrast range of [-80 80]
dotLum = 30;
ringLum = 60; % Max appears to be 94 for EMM302 Display++

% Sizes
% 0.77

%dotSizes = [30 30];
dotSizes = [1.16 1.16].*pxPerDeg(3); % 1.06 degrees
%ringSizes = dotSizes + [10 30];
ringSizes = [1.54 4.72].*pxPerDeg(3);

% Color contrast
nBackgrounds = 17;
maxBackground = 80; % At 45 cd/m2, max green == 88, max red == 103
backgroundContrasts = linspace(-maxBackground, maxBackground, nBackgrounds);
% Remove 0
% backgroundContrasts(backgroundContrasts==0) = [];
% nBackgrounds = nBackgrounds-1;
contrastStep = 1;

% Create trivals
%mwDot = trival('MW', [0 0 dotLum; 0 0 dotLum]);
mwGrayDot = trival('MW', [0 0 dotLum]);
mwRing = trival('MW', [0 0 ringLum; 0 0 ringLum]);
dotContrastRange = 40; % Randomize range for dot

% Positions
viewSpan = 12./2; % Degrees
ringPos = [resolution(1)./2 - viewSpan.*pxPerDeg(3), ...
    resolution(2)./2; ...
    resolution(1)./2 + viewSpan.*pxPerDeg(3), ...
    resolution(2)./2]';
% ringPos = [1920*0.39, 1080*0.5;
%     1920*0.61, 1080*0.5]';

%% Create derived colors
bkgdTri = [backgroundContrasts', ...
    zeros(size(backgroundContrasts))', ...
    ones(size(backgroundContrasts))' .* backgroundLum];
mwBackground = trival('MW', bkgdTri);
[rgbBackground, ff] = mw2rgb(mwBackground);
if any(ff,'all')
    sca
    error('One or more background luminances are oversaturated. Lower luminance.')
end
rgbGrayDot = mw2rgb(mwGrayDot);
[rgbRing, ff_ring] = mw2rgb(mwRing);
if any(ff_ring,'all')
    sca
    error('Ring luminances are oversaturated. Lower luminance.')
end

%% Setup
% Create data table
% t = table([backgroundContrasts, backgroundContrasts, backgroundContrasts, backgroundContrasts]', ...
%     [rgbBackground.Value(:,1); rgbBackground.Value(:,1); rgbBackground.Value(:,1); rgbBackground.Value(:,1)], ...
%     [rgbBackground.Value(:,2); rgbBackground.Value(:,2); rgbBackground.Value(:,2);rgbBackground.Value(:,2)], ...
%     [rgbBackground.Value(:,3); rgbBackground.Value(:,3); rgbBackground.Value(:,3); rgbBackground.Value(:,3)],...
%     [ones(nBackgrounds, 1); ones(nBackgrounds, 1); ones(nBackgrounds,1).*-1; ones(nBackgrounds,1).*-1], ...
%     [ones(nBackgrounds, 1).*1; ones(nBackgrounds, 1).*2; ones(nBackgrounds, 1).*1; ones(nBackgrounds, 1).*2], ...
%     zeros(nBackgrounds*4,1), ...
%     'VariableNames', {'BackContrast','backR','backG','backB','Side','SmallRingIdx','Setting'});
t = table([backgroundContrasts, backgroundContrasts]', ...
    [rgbBackground.Value(:,1); rgbBackground.Value(:,1)], ...
    [rgbBackground.Value(:,2); rgbBackground.Value(:,2)], ...
    [rgbBackground.Value(:,3); rgbBackground.Value(:,3)],...
    [ones(nBackgrounds, 1); ones(nBackgrounds,1).*-1], ...
    [ones(nBackgrounds, 1).*2; ones(nBackgrounds, 1).*1], ...
    zeros(nBackgrounds*2,1), ...
    'VariableNames', {'BackContrast','backR','backG','backB','Side','SmallRingIdx','Setting'});
% Randomize table
t=t(randperm(size(t,1)),:);

instr = sprintf('Adjust the color of the CENTER on the side surrounded by a thinner ring.\nPress any key to start.');
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
setKey = KbName('Space');

% Hide the mouse cursor
% HideCursor;


%% Experiment
[win, rect] = PsychImaging('OpenWindow', screenNum, [0.5 0.5 0.5]);
Screen('BlendFunction', win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
Screen('Flip', win);
ifi = Screen('GetFlipInterval', win);
Screen('TextSize', win, 40);
topPriorityLevel = MaxPriority(win);
if ~debug
    Priority(topPriorityLevel);
    HideCursor;
end

DrawFormattedText(win, instr, 'center', 'center', [1 1 1]);
Screen('Flip', win)
KbWait;
WaitSecs(0.5);

% Experiment block
for i = 1:size(t,1)
    % Get current trial
    currTrial = t(i,:);
    % Randomize starting dot contrast at start of each trial
    randStartContrast = round(dotContrastRange - (rand*dotContrastRange*2));
    currContrast = randStartContrast;
    mwCurrDot = trival('MW', [currContrast, 0 , dotLum]);

    % Audio cue for trial
    if currTrial.Side == -1
        Beeper(2000);
    elseif currTrial.Side == 1
        Beeper(6000);
    end

    % Assign ring to correct side
    if currTrial.SmallRingIdx == 2
        currRingSizes = fliplr(ringSizes);
    else
        currRingSizes = ringSizes;
    end

    % Start trial
    trialEnded = 0;
    while ~trialEnded
        % Reset dot contrast, generate RGB
        mwCurrDot.Value(1) = currContrast;
        rgbDot = mw2rgb(mwCurrDot);

        % Assign dot to correct side
        if currTrial.Side == -1
            leftDot = rgbDot';
            rightDot = rgbGrayDot';
        elseif currTrial.Side == 1
            leftDot = rgbGrayDot';
            rightDot = rgbDot';
        end

        % Draw background
        rgbBackground = [t{i,"backR"}, t{i,"backG"}, t{i,"backB"}];
        Screen('FillRect', win, rgbBackground);

        % Draw rings
        Screen('DrawDots', win, ringPos, currRingSizes, rgbRing.Value', [], 1);
        Screen('DrawDots', win, ringPos(:,1), dotSizes(1), leftDot.Value', [], 1);
        Screen('DrawDots', win, ringPos(:,2), dotSizes(2), rightDot.Value', [], 1);
        if debug
            DrawFormattedText(win, sprintf('Adjust %d\nBkgd = %.4f\nContrast = %.4f', ...
                currTrial.Side,currTrial.BackContrast, currContrast), 'center', rect(4)/4, [1 1 1]);
        end
        Screen('Flip', win);

        % Get response

        [keyIsDown,secs, keyCode] = KbCheck(-1);
        if keyCode(escapeKey)
            cd(pwd);
            ShowCursor;
            sca;
            return
        elseif keyCode(leftKey)
            %beep
            currContrast = currContrast - contrastStep;
        elseif keyCode(rightKey)
            %beep
            currContrast = currContrast + contrastStep;
        elseif keyCode(setKey)
            Beeper(600)
            t.Setting(i) = currContrast;
            trialEnded = 1;
        end
        if currContrast > maxBackground
            %Beeper(8000)
            beep
            currContrast = maxBackground;
        elseif currContrast < -maxBackground
            %Beeper(8000)
            beep
            currContrast = -maxBackground;
        end
        % Wait for end of response
        %WaitSecs(0.1);

        %    endResp = GetSecs;
        %    rt = endResp - startResp;
    end
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
writetable(sortrows(t(:, ["BackContrast", "Side", "SmallRingIdx", "Setting"])), subjID + ".xlsx", 'Range', "A1:D" + string(size(t,1)+1));
writetable(tAvgSmall, subjID + ".xlsx", 'Range', "F1:G" + size(tAvgSmall,1)+1);
writetable(tSDSmall(:,2), subjID + ".xlsx", 'Range', "H1:H" + size(tSDSmall,1)+1);
writetable(tAvgLarge(:,2), subjID + ".xlsx", 'Range', "I1:I" + size(tAvgLarge,1)+1);
writetable(tSDLarge(:,2), subjID + ".xlsx", 'Range', "J1:J" + size(tSDLarge,1)+1);

% tOut = [tAdjSmall, tAdjLarge(:,2)];
% writetable(tOut, subjID + ".xlsx", 'Range', "A1:C" + string(size(tOut,1)+1));
% writecell({'Mean', 'SD'}, subjID + ".xlsx", 'Range', "E1:F1");
% writematrix([avg, sd], subjID + ".xlsx", 'Range', "E2:F"+ string(size(avg,1)+1));

ShowCursor;
sca
