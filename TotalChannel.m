classdef TotalChannel
    
    properties (GetAccess='public', SetAccess='public')
        SymbolRate           % time precision of analog wave
        SampleRate              % AD/DA sample rate
        
        ConvCode                % Convolutional Coding
        
               
        
    end
    
    methods
        function obj = TotalChannel(SymbolRate1,SampleRate1,ConvCode1)
            obj.SymbolRate = 2000;
            obj.SampleRate = 32000;
            obj.ConvCode = [];
            if nargin >= 1                
                obj.SymbolRate = SymbolRate1;
                if nargin >= 2
                    obj.SampleRate = SampleRate1;
                    if nargin >= 3
                        obj.ConvCode = ConvCode1;
                    end
                end
            end
            
        end
        
        function [ErrorRate,EbByN0,SPDgram] = RunWithoutEncoding(obj,DataIn,N0)
            %N0 sigle side SPD of noise
            % 4QAM
            
            fs = obj.SampleRate;
            rs = obj.SymbolRate;
            alpha = 0.5;
            DataIn = reshape(DataIn,[],1);
            Vol = reshape(DataIn,2,[])';
            Vol = bi2de(Vol);
            Vol = bin2gray(Vol,'psk',4);
            Vol = exp(1j*(2*Vol)/4*pi);
            Tsymbol = (0:length(Vol)-1)'/rs;            
            T0 = (length(Vol)+200)/rs;
            
            Tg = 50/rs;     % G(t) duration
            Ng = fs*Tg;         
            freq = [0:ceil(Ng/2)-1,-floor(Ng/2):-1]'/Tg;
            
            G = sqrt((abs(freq)<=(1-alpha)*rs/2)+(abs(freq)/rs*2>1-alpha&abs(freq)/rs*2<1+alpha).*cos((abs(freq)/rs*2-1+alpha)*2*pi/8/alpha).^2);
            g = real(ifft(G));
%             timeg = [0:ceil(Ng/2)-1 -floor(Ng/2):-1]'/fs;
            g = circshift(g,round(Tg/2*fs))/sqrt(sum(abs(g).^2));
            
            timew = (1:T0*fs)'/fs;
%             waveBBT = pulstran(timew,[Tsymbol,Vol],g,fs);   %transmitted wave on baseband

            waveBBT = zeros(T0*fs,1);
            waveBBT(((1:length(Vol))-1)*fs/rs+1) = Vol; % pulse train
            waveBBT = conv(waveBBT,g,'full');
            waveBBT(T0*fs+1:end) = [];
            plotc(waveBBT(1:fs/rs*20))
            
            SPDgram = [ [0:ceil(T0*fs/2)-1 -floor(T0*fs/2):-1]' abs(fft(waveBBT)).^2];
            
            waveT = real(waveBBT.*exp(2j*pi*1850*timew));
            
            noise = random('Normal',0,N0,size(waveT));
            
            waveR = waveT+noise;
            waveBBR = 2*waveR.*exp(-2j*pi*1850*timew);
            g1 = ifft(fft(g)');
            waveMF = conv(waveBBR,g1,'full');
            
            plotc(waveMF(1:30*fs/rs))
            
            VolR = waveMF(length(g)+1:fs/rs:length(waveMF));
            Vol0 = bin2gray(0:3,'psk',4);
            Vol0 = exp(1j*(2*Vol0)/4*pi);
            est = -abs(VolR-Vol0).^2;
            [~,symbolR] = max(est,[],2);
            symbolR = symbolR-1;
            symbolR = de2bi(symbolR);
            symbolR = reshape(symbolR',[],1);
            symbolR = symbolR(1:length(DataIn));
            ErrorRate = mean(abs(symbolR-DataIn));
        end
    end
    

end
    
    
