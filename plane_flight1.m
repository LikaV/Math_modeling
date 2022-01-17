 
    DT = 0.1;
    
    u = 15.041/60/180*pi;
    e2 = 6.694e-3;
    g = 9.8;
    Vot = 21.54;
    m = 1912;
    p = 97e4;
    Hb =10.7;
    Vb = 26.6;
    Xa=21.54*22.44*1.2*1.5;
    Lp = 1000*Vot^2/(2*g*(p/(m*g) - 0.2));
    Lph = 1000*(m*g/(p - Xa/10))*((Vb^2 - Vot^2)/(2*g) + Hb);
    T = [0 : 0.1 : Lp/Vot]';
    Ge = 9.7803253359;
    Gp = 9.8321849378;
    
    gamma = T*0;
    theta = T*0;
    psi0 = atan(17.5/16.5);
    psi = psi0+0*T;
    phi0 = (11+32/60+13.47/3600)/180*pi;
    lam0 = (104+50/60+03.58/3600)/180*pi;

    rn = 6378137*(1-e2)/(1 - e2*sin(phi0))^(3/2);
    re = 6378137/(1 - e2*sin(phi0))^(1/2);
    H0 = 13 + 3; %высота относительно уровня моря
    H = H0 + T*0; %высота в каждый момент времени матрица размера как и Т ,но с числами H0
    x = T/Lp*Vot;
    
    V0 = 0;
    Ver = ((V0 - Vot)/2*cos(pi*x) + (V0 + Vot)/2)*sin(psi0);
    Vnr = ((V0 - Vot)/2*cos(pi*x) + (V0 + Vot)/2)*cos(psi0);
    V = [Ver, Vnr, T*0];
    Vdifer = -pi*(V0 - Vot)/2*sin(pi*x)*sin(psi0);
    Vdifnr = -pi*(V0 - Vot)/2*sin(pi*x)*cos(psi0);
    Vdif = [Vdifer, Vdifnr, T*0];
    Sgamma = [1, 0, 0; 0, 0, 1; 0, -1, 0];
    Spsi = [sin(psi0), cos(psi0), 0; -cos(psi0), sin(psi0), 0; 0, 0, 1];
    Sr = Sgamma*Spsi;
    phi(1, 1) = phi0;
    lam(1, 1) = lam0;
    G(1, :) = [0, 0, -(Ge*(cos(phi0))^2 + Gp*(1 - e2)^(1/2)*(sin(phi0))^2)/(1-e2*(sin(phi0))^2)^(1/2)];
    Fz(:, 1) = Sr*(- G(1, :)');
    for i = 2 : size(T, 1)
        phi(i, 1) = phi(i - 1, 1) + Ver(i, 1)*DT/(rn + H0);
        lam(i, 1) = lam(i - 1, 1) + Vnr(i, 1)*DT/((re + H0)*cos(phi(i, 1)));
        G(i, :) = [0, 0, -(Ge*(cos(phi(i, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phi(i, 1)))^2)/(1-e2*(sin(phi(i, 1)))^2)^(1/2)];
        Omega = [-V(i, 2)/(rn+H0), V(i, 1)/(re+H0), V(i, 1)*tan(phi(i, 1))/(re+H0)];
        Ux = [0, u*cos(phi(i, 1)), u*sin(phi(i, 1))];
        Omega_ = [0, Omega(3), -Omega(2); Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
        Ux_ = [0, Ux(3), -Ux(2); Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
        Fz(:, i) = Sr*(Vdif(i, :)' - (Omega_ + 2*Ux_)*V(i, :)' - G(i, :)');
    end
    FZr = [T, Fz']
   % save fz.txt FZ -ascii
    
    phi1 = phi(end, 1);
    lam1 = lam(end, 1);

    flight0 = [T lam*180/pi phi*180/pi H gamma theta psi];
  
    
    T1 = [Lp/Vot + 0.1 : 0.1 : Lp/Vot + 5]';  
    x = (T1 - Lp/Vot)/5;
    H1 = H0 + Hb;
    gammaot = T1*0;
    psiot = psi0+0*T1;
    v0 = 0;
    v1 = 2;
   
    theta1 = 15/180*pi;
    Veot = ((Vot - Vb)/2*cos(pi*x) + (Vb + Vot)/2)*sin(psi0);
    Vnot = ((Vot - Vb)/2*cos(pi*x) + (Vb + Vot)/2)*cos(psi0);
    Vup = (v0 - v1)/2*cos(pi*(1-x)) + (v0 + v1)/2;
    V2 = [Veot, Vnot, Vup];
    V2difer = -pi*(Vot - Vb)/2*sin(pi*x)*sin(psi0);
    V2difnr = -pi*(Vot - Vb)/2*sin(pi*x)*cos(psi0);
    V2difup = pi*(v0 - v1)/2*sin(pi*(1 - x));
    V2dif = [V2difer, V2difnr, V2difup];
    
    thetaot = - theta1/2*cos(pi*x) + theta1/2;
    phiot(1, 1) = phi1;
    lamot(1, 1) = lam1;
    Hot(1, 1) = H0;
    Got(1, :) = [0, 0, -(Ge*(cos(phiot(1, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiot(1, 1)))^2)/(1-e2*(sin(phiot(1, 1)))^2)^(1/2)];
    Omega = [-V2(1, 2)/(rn + H0), V2(1, 1)/(re + H0), V2(1, 1)*tan(phiot(1, 1))/(re + Hot(1, 1))];
    Ux = [0, u*cos(phiot(1, 1)), u*sin(phiot(1, 1))];
    Omega_ = [0, Omega(3), -Omega(2); Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
    Ux_ = [0, Ux(3), -Ux(2); Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
    Fzot(:, 1) = Sr*(V2dif(1, :)' - (Omega_ + 2*Ux_)*V2(1, :)' - Got(1, :)');
    for i = 2 : size(T1, 1)
        Hot(i, 1) = Hot(i - 1, 1) + Veot(i, 1)*DT;
        phiot(i, 1) = phiot(i - 1, 1) + Vnot(i, 1)*DT/(rn + Hot(i, 1));
        lamot(i, 1) = lamot(i - 1, 1) + Vup(i, 1)*DT/((re + Hot(i, 1))*cos(phiot(i, 1)));
        Got(i, :) = [0, 0, -(Ge*(cos(phiot(i, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiot(i, 1)))^2)/(1-e2*(sin(phiot(i, 1)))^2)^(1/2)];
        Omega = [-V2(i, 2)/(rn + Hot(i, 1)), V2(i, 1)/(re + Hot(i, 1)), V2(i, 1)*tan(phiot(i, 1))/(re + Hot(i, 1))];
        Ux = [0, u*cos(phiot(i, 1)), u*sin(phiot(i, 1))];
        Omega_ = [0, Omega(3), -Omega(2); Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
        Ux_ = [0, Ux(3), -Ux(2); Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
        Stheta = [cos(thetaot(i, 1)), 0, sin(thetaot(i, 1)); 0, 1, 0; -sin(thetaot(i, 1)), 0, cos(thetaot(i, 1))];
        S = Sgamma*Stheta*Spsi;
        Fzot(:, i) = S*(V2dif(i, :)' - (Omega_ + 2*Ux_)*V2(i, :)' - Got(i, :)');
    end
    FZot = [T1, Fzot'];
    %FZ = [FZr; FZot];
    
    flight1 = [T1 lamot*180/pi phiot*180/pi Hot gammaot thetaot psiot];
    
    Rn = rn;
    Re = re;
    dt = 0.1;
    T2 = [T1(end,1):0.1:T1(end,1)+180]';
    Vkr = 296/3600*1000; 
    x = (T2-T1(end,1))/180; 
    Veup = ((Vb - Vkr)/2*cos(pi*x) + (Vkr + Vb)/2)*sin(psi0); 
    Vnup = ((Vb - Vkr)/2*cos(pi*x) + (Vkr + Vb)/2)*cos(psi0); 
    Vup = (v0 - v1)/2*cos(pi*x) + (v0 + v1)/2; 
    V3 = [Veup, Vnup, Vup];
    V3difer = -pi*(Vb - Vkr)/2*sin(pi*x)*sin(psi0);
    V3difnr = -pi*(Vb - Vkr)/2*sin(pi*x)*cos(psi0);
    V3difup = -pi*(v0 - v1)/2*sin(pi*x);
    V3dif = [V3difer, V3difnr, V3difup];
    %gammaup = (0 - 30/180*pi)/2*cos(pi*x.^10) + (30/180*pi + 0)/2; 
    gammaup = 0*T2;
    thetaup = thetaot(end, 1)/2*cos(pi*x) + thetaot(end, 1)/2; 
    psiup = psi0+0*T2; 
    phiup(1, 1) = phiot(end, 1); 
    lamup(1, 1) = lamot(end, 1); 
    Hup(1, 1) = Hot(end, 1); 
    Gup(1, :) = [0, 0, -(Ge*(cos(phiot(1, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiot(1, 1)))^2)/(1-e2*(sin(phiot(1, 1)))^2)^(1/2)];
    Omega = [-V3(1, 2)/(rn + Hup(1, 1)), V3(1, 1)/(re + Hup(1, 1)), V3(1, 1)*tan(phiot(1, 1))/(re + Hup(1, 1))];
    Ux = [0, u*cos(phiup(1, 1)), u*sin(phiup(1, 1))];
    Omega_ = [0, Omega(3), -Omega(2); Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
    Ux_ = [0, Ux(3), -Ux(2); Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
    Stheta = [cos(thetaup(1, 1)), 0, sin(thetaup(1, 1)); 0, 1, 0; -sin(thetaup(1, 1)), 0, cos(thetaup(1, 1))];
    S = Sgamma*Stheta*Spsi;
    Fzup(:, 1) = S*(V3dif(1, :)' - (Omega_ + 2*Ux_)*V3(1, :)' - Gup(1, :)');
    for i = 2 : size(T2, 1) 
         Hup(i, 1) = Hup(i - 1, 1) + Vup(i, 1)*dt; 
         phiup(i, 1) = phiup(i - 1, 1) + Vnup(i, 1)*dt/(Rn + Hup(i, 1)); 
         lamup(i, 1) = lamup(i - 1, 1) + Veup(i, 1)*dt/((Re + Hup(i, 1))*cos(phiup(i, 1)));
         Gup(i, :) = [0, 0, -(Ge*(cos(phiup(i, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiup(i, 1)))^2)/(1-e2*(sin(phiup(i, 1)))^2)^(1/2)];
        Omega = [-V3(i, 2)/(rn + Hup(i, 1)), V3(i, 1)/(re + Hup(i, 1)), V3(i, 1)*tan(phiup(i, 1))/(re + Hup(i, 1))];
        Ux = [0, u*cos(phiup(i, 1)), u*sin(phiup(i, 1))];
        Omega_ = [0, Omega(3), -Omega(2); Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
        Ux_ = [0, Ux(3), -Ux(2); Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
        Stheta = [cos(thetaup(i, 1)), 0, sin(thetaup(i, 1)); 0, 1, 0; -sin(thetaup(i, 1)), 0, cos(thetaup(i, 1))];
        S = Sgamma*Stheta*Spsi;
        Fzup(:, i) = S*(V3dif(i, :)' - (Omega_ + 2*Ux_)*V3(i, :)' - Gup(i, :)');
    end 
    FZup = [T2, Fzup'];
    FZ = [FZr; FZot; FZup];
    save fz.txt FZ -ascii
    flight2 = [T2 lamup*180/pi phiup*180/pi Hup gammaup thetaup psiup]; 
    flight = [flight0 ; flight1; flight2 ];
    save flight.txt flight -ascii
    tgeo2kml(flight)
    
              