function omega = objective(x)
    rho = 997;
    I = (1/12)*0.85*0.5*0.5;
    
    omega(1) = 0;
    t(1) = 0;
    alpha(1) = 0;
    f_total(1) = 0;
    j=1;
    Kn = 1;
    while sum(omega)==0 || round(omega(j),5) ~= round(omega(j-1),5)
        total = 0;
        for i = 1:x(2)
            force = Kn*(rho*pi/(2*x(2)))*(((x(6))*x(3))^2-omega(j)*x(6)*x(3)^3);
            distance = (x(4)*(i-1)+x(5));
            total = (force-omega(j)*0.05)*distance;
        end
        f_total(j+1) = total/x(2);
        alpha(j+1) = (2*total-omega(j)*0.01)*I;
        omega(j+1) = omega(j)+alpha(j+1)*t(j);
        t(j+1)= t(j)+1;
        j=j+1;
    end
    omega = -omega(j);
    
