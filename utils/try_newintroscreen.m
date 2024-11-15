%% Preparation scripts

% Screen Setup
PsychDefaultSetup(2);
% Get the screen numbers
screens = Screen('Screens');
% Select the external screen if it is present, else revert to the native
% screen
screenNumber = max(screens);
% Define black, white and grey
black = BlackIndex(screenNumber);
white = WhiteIndex(screenNumber);
grey = white / 2;
% Open an on screen window and color it grey
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);

% Set the blend funnction for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% Get the size of the on screen window in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
ifi = Screen('GetFlipInterval', window);
% Get the centre coordinate of the window in pixels
[xCenter, yCenter] = RectCenter(windowRect);
% Set the text size
Screen('TextSize', window, 50);


%% Experimental scripts
cue = '1 2';
topText = 'You will hear two sounds (1 and 2). Keep them in mind. \n Then you will see a cue on screen: ';
bottomText = 'After a delay, you will see a green circle. \n Repeat the two sounds as the green circle appears.';

Screen('TextSize', window, 100);
DrawFormattedText(window, cue, 'center', 'center', [1 1 1]);

[normBoundsRect, ~] = Screen('TextBounds', window, cue);
cueHeight = normBoundsRect(4) - normBoundsRect(2);

Screen('TextSize', window, 50);
DrawFormattedText(window, topText, 'center', yCenter - cueHeight / 2 - 120, [1 1 1]);

DrawFormattedText(window, bottomText, 'center', yCenter + cueHeight / 2 + 120, [1 1 1]);

Screen('Flip', window);

WaitSecs(10);

sca;