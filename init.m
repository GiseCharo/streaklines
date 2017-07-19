clear;close all;clc;
dbstop if error;

tests=[ pwd '/test' ];
fct = genpath([ pwd '/functions' ]);
frt = genpath([ pwd '/main' ]);
addpath(fct)
addpath(frt)
addpath(tests);
clear current_tests mains
home;
