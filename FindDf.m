function df = FindDf(ConvCode)
    nm = ConvCode.nm;
    m = nm(2)-1;
    df = Inf;
    for nn = 1:2^(2*m)
        code = ConvCode.ConvEncoder([de2bi(nn)]');
        df = min(sum(code),df);
        
    end

end