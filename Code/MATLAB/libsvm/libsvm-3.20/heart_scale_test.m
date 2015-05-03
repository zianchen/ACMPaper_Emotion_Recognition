%% HowToClassifyUsingLibsvm
% by faruto @ faruto's Studio~
% http://blog.sina.com.cn/faruto
% Email:faruto@163.com
% http://www.matlabsky.com
% http://www.mfun.la
% http://video.ourmatlab.com
% last modified by 2010.12.27
%% a litte clean work
tic;
close all;
clear;
clc;
format compact;
%% 

% ??????
[heart_scale_label,heart_scale_inst] =libsvmread('heart_scale'); 
data = heart_scale_inst;
label = heart_scale_label;

% ???200???????????70?????????
ind = 200;
traindata = data(1:ind,:);
trainlabel = label(1:ind,:);
testdata = data(ind+1:end,:);
testlabel = label(ind+1:end,:);

% ????????????
model = svmtrain(trainlabel,traindata,'-s 0 -t 2 -c 1.2 -g 2.8');

% ????model??
model
Parameters = model.Parameters
Label = model.Label
nr_class = model.nr_class
totalSV = model.totalSV
nSV = model.nSV 

% 
[ptrain] = svmpredict(trainlabel,traindata,model);

% 
[ptest] = svmpredict(testlabel,testdata,model);

%%
toc;