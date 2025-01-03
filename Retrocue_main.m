% Main script to run the Retrocue experiment:

function Retrocue_main(subject,practice,startblock)
% subject (string) = subject number, string
% practice (integer) =  practice or full task. 
%                  1 = full practice + real task (starting from block 1).
%                  2 = mixed session in practice only + real task (starting from block 1).
%                  other number = real task (starting from startblock).
% startblock (integer) = starting block at the beginning.

cd utils

if practice == 1
    if startblock~=1
        disp('In a mode containing Practice session, startblock can only be 1. Changed to 1 automatically');
    end
    run_code=Retro_Repeat_WithinBlock_util(subject,1,1);
    while run_code~=-1
       run_code=Retro_Repeat_WithinBlock_util(subject,run_code,1);
    end

elseif practice == 2
    if startblock~=1
        disp('In a mode containing Practice session, startblock can only be 1. Changed to 1 automatically');
    end
    run_code=Retro_Repeat_WithinBlock_util(subject,2,1);
    while run_code~=-1
      run_code=Retro_Repeat_WithinBlock_util(subject,run_code,1);
    end
    
else

    Retro_Repeat_WithinBlock_util(subject,0,startblock);
end

cd ..