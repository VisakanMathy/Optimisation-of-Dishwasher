function out = checkForInt(x)
    if((x - round(x))^2 < 0.1)
        out = 0;
    else
        out = 1;
    end