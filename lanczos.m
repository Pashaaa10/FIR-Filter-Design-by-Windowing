function w = lanczos(M,L)
    n = 0:M-1;
    w = (sin(2*pi*((2*n-M+1)/(2*(M-1))))./(2*pi*((2*n-M+1)/(2*(M-1))))).^L;
    w = w'; %transpunere
end