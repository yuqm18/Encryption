clear;
close all


Channel = TotalChannel();
Channel.ConvCode = ConvCode([11 13 15]);


%data = repmat([1;zeros(1000,1)],10,1);
nN0 = 20;
rep = 50;
EbByN0 = 10.^linspace(-1,1,nN0);
ErrorRate = zeros(nN0,rep);
for n = 1:rep
    data = rand(8192,1)<0.5;
    for m = 1:nN0
        [ErrorRate(m,n)] = Channel.Run(data,EbByN0(m));
    end
end
ErrorRateLine = mean(ErrorRate,2);
EbByN0dB = 10*log10(EbByN0);
semilogy(EbByN0dB,ErrorRateLine,'*')
hold on
semilogy(EbByN0dB,normcdf(-sqrt(2*EbByN0)))
hold off
title('未编码')
ylabel('ErrorRate')
xlabel('Eb/N0 \\dB')
legend('仿真点','理论值')



