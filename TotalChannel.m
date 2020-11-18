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
            
            Vol = reshape(DataIn,[],2);
            Vol = bi2de(Vol);
            Vol = bin2gray(Vol,'psk',4);
            Vol = exp(1j*(1+2*Vol)/4*pi);
            Tsymbol = (0:length(Vol)-1)'/rs;            
            T0 = (length(Vol)+20)/rs;
            
            Tg = 20/rs;     % G(t) duration
            Ng = fs*Tg;         
            freq = [0:ceil(Ng/2)-1,-floor(Ng/2):-1]'/Tg;
            
            G = sqrt((abs(freq)<=(1-alpha)*rs)+(abs(freq)/rs>1-alpha&abs(freq)/rs<1+alpha).*cos((abs(freq)/rs-1+alpha)*2*pi/8/alpha).^2);
            g = real(ifft(G));
%             timeg = [0:ceil(Ng/2)-1 -floor(Ng/2):-1]'/fs;
            g = circshift(g,round(Tg/2*fs))/rms(g);
            
            timew = (1:T0*fs)'/fs;
            waveBBT = pulstran(timew,[Tsymbol,Vol],g,fs);   %transmitted wave on baseband
            plotc(waveBBT(1:fs/rs*20))
            
            SPDgram = [ [0:ceil(T0*fs/2)-1 -floor(T0*fs/2):-1]' abs(fft(waveBBT)).^2];
            
            waveT = real(waveBBT.*exp(2j*pi*1850*timew));
            
            noise = random('Normal',0,N0,size(waveT));
            
            waveR = waveT+noise;
            waveBBR = 2*waveR.*exp(-2j*pi*1850*timew);
            g1 = ifft(fft(g)');
            waveMF = conv(waveBBR,g1,'valid');
            
            plotc(waveMF(1:10*fs/rs))
            
            
        end
    end
    

end
    
    
