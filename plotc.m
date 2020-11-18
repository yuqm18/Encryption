function plotc(varargin)
if nargin == 1
    y = varargin{1};
reshape(y,[],1);
plot([real(y) imag(y)]);
else
    y = varargin{2};
    x = varargin{1};
    reshape(y,[],1);
plot(x,[real(y) imag(y)]);
end

end

