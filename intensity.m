%%% I = A*exp(beta*T) %%%
function I = intensity(T)
    I0 = 55; %dB
    IF = 70;
    TIME = 3.1;
    beta = log(IF/I0)/TIME;
    I = I0*exp(beta*T);
end