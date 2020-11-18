clear;
close all

data = rand(8192,1)<0.5;
Channel = TotalChannel();

Channel.RunWithoutEncoding(data,1);