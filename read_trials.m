function [trial_No,block_No,Syll1_No,Syll2_No,Retro_No,Retro_Brightness]=read_trials(subject,stim_Tags,retro_Tags)

filename = fullfile('trials',[subject '_trial_list.xlsx']);

if ~isfile(filename)
    backupFiles = dir(fullfile('trials', '*backup_trial_list*.xlsx'));
    
    if ~isempty(backupFiles)
        backupFile = fullfile('trials', backupFiles(1).name);
        movefile(backupFile, filename);
        disp(['File renamed from ' backupFiles(1).name ' to ' subject '_trial_list.xlsx']);
    else
        error('No backup file found in the directory.');
    end

end

data = readtable(filename);
disp('File loaded successfully.');


data = readtable(filename);
numRows = height(data);

trial_No =  nan(numRows, 1);
block_No =  trial_No;
Syll1_No =  trial_No;
Syll2_No =  trial_No;
Retro_No =  trial_No;
Retro_Brightness =  trial_No;

for i = 1:numRows

    trial_No(i) = data.Trial(i);
    block_No(i) = data.Suggested_Block(i);
    Retro_Brightness(i) = 0.2+0.8*data.Cue_brightness(i); % Brightness ranging from 0.2 to 1

    Syll1 = data.Syllable_1(i);
    switch Syll1{1}
        case stim_Tags{1}
            % /i/
            Syll1_No(i) =  1;
        case stim_Tags{2}
            % /o/
            Syll1_No(i) =  2;
        case stim_Tags{3}
            % /a/
            Syll1_No(i) =  3;
    end

    Syll2 = data.Syllable_2(i);
    switch Syll2{1}
        case stim_Tags{1}
            % /i/
            Syll2_No(i) =  1;
        case stim_Tags{2}
            % /o/
            Syll2_No(i) =  2;
        case stim_Tags{3}
            % /a/
            Syll2_No(i) =  3;
    end

    Retro =  data.Retrocue(i);
    switch Retro{1}
        case retro_Tags{1}
            % "REP_BTH"
            Retro_No(i) =  1;
        case retro_Tags{2}
            % "REV_BTH"
            Retro_No(i) =  2;
        case retro_Tags{3}
            % "REP_1ST"
            Retro_No(i) =  3;
        case retro_Tags{4}
            % "REP_2ND"
            Retro_No(i) =  4;
        case retro_Tags{5}
            % "DRP_BTH"
            Retro_No(i) =  5;
    end

end