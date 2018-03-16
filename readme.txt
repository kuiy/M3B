
The program was successfully compiled in VC++ 6.0 and run at command window at Windows 7.


m3b alarm20000.data 20000 37 25 0.01


where

alarm.data: trainding data including 20000 data samples (rows) and 37 variables (columns)

20000 is the number of data cases in alarm.data,

37 is the number of variables in alarm.data,

25: the index of the target variable in alarm.data (in C++, since the index of an array is from 0, 25 represents the 26th variables in 
alarm.data),

0.01 is the significance level.




output:

result_mb.txt: store the mb discovered

result_pc.txt: store the pc discovered