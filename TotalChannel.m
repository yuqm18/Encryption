classdef TotalChannel
    
    properties (GetAccess='public', SetAccess='public')
        SymbolRate           % time precision of analog wave
        SampleRate              % AD/DA sample rate
        
        ConvCode                % Convolutional Coding
        
        ModType                 % Modulation Type 'psk','qam'
               
        
    end
    
    methods
        function obj = TotalChannel(SymbolRate1,SampleRate1,ConvCode1)
            obj.SymbolRate = 2000;
            obj.SampleRate = 32000;
            obj.ConvCode = [];
            obj.ModType = [];
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
        
        function [ErrorRate,SPDgram] = Run(obj,DataIn,EbByN0)
            %N0 sigle side SPD of noise
            
            
            fs = obj.SampleRate;
            rs = obj.SymbolRate;
            alpha = 0.5;

            if ~isempty(obj.ConvCode)
                DataIn1 = obj.ConvCode.ConvEncoder(DataIn);
                M = obj.ConvCode.nm(1);
            else 
                DataIn1 = DataIn;
                M = [];
            end
            
            [Vol,MVol,StdVol] = obj.LevelModulation(DataIn1,M);
            Tsymbol = (0:length(Vol)-1)'/rs;            
            T0 = (length(Vol)+200)/rs;
            
            Tg = 50/rs;     % G(t) duration
            Ng = fs*Tg;         
            freq = [0:ceil(Ng/2)-1,-floor(Ng/2):-1]'/Tg;
            
            G = sqrt((abs(freq)<=(1-alpha)*rs/2)+(abs(freq)/rs*2>1-alpha&abs(freq)/rs*2<1+alpha).*cos((abs(freq)/rs*2-1+alpha)*2*pi/8/alpha).^2);
            g = real(ifft(G));
%             timeg = [0:ceil(Ng/2)-1 -floor(Ng/2):-1]'/fs;
            g = circshift(g,round(Tg/2*fs))/sqrt(sum(abs(g).^2));
            Eb = sum(abs(g).^2)/fs...   % energy of g(t)
                *mean(abs(StdVol).^2)...    % average Level，equal prob.
                *length(Vol)/length(DataIn)/2;
            
            timew = (1:T0*fs)'/fs;
%             waveBBT = pulstran(timew,[Tsymbol,Vol],g,fs);   %transmitted wave on baseband

            waveBBT = zeros(T0*fs,1);
            waveBBT(((1:length(Vol))-1)*fs/rs+1) = Vol; % pulse train
            waveBBT = conv(waveBBT,g,'full');
            waveBBT(T0*fs+1:end) = [];
            plotc(waveBBT(1:fs/rs*20))
            
            SPDgram = [ [0:ceil(T0*fs/2)-1 -floor(T0*fs/2):-1]' abs(fft(waveBBT)).^2];
            
            waveT = real(waveBBT.*exp(2j*pi*1850*timew));
            
            noise = random('Normal',0,sqrt(Eb/EbByN0*fs/2),size(waveT));
            
            waveR = waveT+noise;
            waveBBR = 2*waveR.*exp(-2j*pi*1850*timew);
            g1 = ifft(fft(g)');
            waveMF = conv(waveBBR,g1,'full');
            
            plotc(waveMF(1:30*fs/rs))
            
            VolR = waveMF(length(g)+1:fs/rs:length(waveMF));
%             Vol0 = bin2gray(0:3,'psk',4);
%             Vol0 = exp(1j*(2*Vol0)/4*pi);
            est = -abs(VolR-StdVol).^2;
            if isempty(obj.ConvCode)
            [~,symbolR] = max(est,[],2);
            symbolR = symbolR-1;
            symbolR = de2bi(symbolR);
            symbolR = reshape(symbolR',[],1);            
            else
                symbolR = obj.ConvCode.ConvDecoder(est);
            end
            symbolR = symbolR(1:length(DataIn));
            ErrorRate = mean(abs(symbolR-DataIn));
        end
        
        function [Vol, M,StdVol] = LevelModulation(obj,DataIn,varargin)
            if nargin<3 || isempty(varargin{1})
                M = max(ceil(log2(length(DataIn)/(5*obj.SymbolRate))),0)+1;
            else
                M = varargin{1};
            end
            if strcmp(obj.ModType,'psk')
                modtype = @pskmod;
            elseif strcmp(obj.ModType,'qam')
                modtype = @qammod;
            else
                modtype = @pskmod;
            end
            DataIn = [reshape(DataIn,[],1);zeros(mod(M-mod(length(DataIn),M),M),1)];
            Vol = reshape(DataIn,M,[])';
            Vol = bi2de(Vol);
            Vol = modtype(Vol,2^M);
            StdVol = modtype(0:2^M-1,2^M);
            
        end
        
    end
    

end
    
    
