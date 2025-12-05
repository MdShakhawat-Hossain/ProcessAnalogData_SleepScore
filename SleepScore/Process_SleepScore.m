clc; clear all; close all;
delete(findall(0,'Type','Figure')); %Close every figure regardless of type

opts = struct('BinSizeSeconds',5,'TrialFolder',pwd);
SleepScore_KECK(opts);