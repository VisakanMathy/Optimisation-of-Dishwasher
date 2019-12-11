function out = checkForInt(x)
    if((x - round(x))^2 > 0.1)
        out = 1;
    else
        out = 0;
    end