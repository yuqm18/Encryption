clear;
close all

data = rand(8192,1)<0.5;
Channel = TotalChannel();

%data = repmat([1;zeros(1000,1)],10,1);


[ErrorRate] = Channel.RunWithoutEncoding(data,1);