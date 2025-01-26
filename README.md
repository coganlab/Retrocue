# Task scripts for Retrocue Experiment

To run this task script, use Matlab:  
`Retrocue_main(subject,practice,startblock)`  

**Modes:**  
`Retrocue_main(subject,1,1)` full practice + real task (starting from block 1).  
`Retrocue_main(subject,2,1)` mixed session in practice only + real task (starting from block 1).  
`Retrocue_main(subject,0,1)` starting block at the beginning.  

**Important:** During the task, in each trial, <font color="red"> patient should repeat continuously **without leaving a gap** between the two syllables </font>. 

![RetroCue task pipeline](stim/key%20press%20instruct.png)      

*To generate random trial lists*  
Clone https://github.com/Baishenliang/lexical_retro_delay_expdesign, install the environment, and fully run the `bsliang_retrocue_expdesign.ipynb`.   
A `backup_trial_list_xxx.xlsx` will be generated and simultaneously sent to the `triallists` folder, and the `trials` folder in the Retrocue_taskscripts local repository (make sure to edit the path beforehead).  
A report of randomization estimation `backup_trial_list_syllablescores_xxx.xlsx` will also be generated and sent to the `triallists` folder.   

Developed by **Baishen Liang** (liangbs95@gmail.com;baishen.liang@duke.edu), adopting the lexical delayed repeat tasks in Cogan Lab.  