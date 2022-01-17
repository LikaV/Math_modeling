 
    DT = 0.1;
    u = 15.041/3600/180*pi;
    e2 = 6.694e-3;
    g = 9.8;
    V0 = 0;
    Vot = 21.54;
    m = 1912;
    p = 97e4;
    Hb =10.7;
    Vb = 26.6;
    Xa=21.54*22.44*1.2*1.5;
    Lp = 1000*Vot^2/(2*g*(p/(m*g) - 0.2));
    Lph = 1000*(m*g/(p - Xa/10))*((Vb^2 - Vot^2)/(2*g) + Hb);
    Ge = 9.7803253359;
    Gp = 9.8321849378;
    H0 = 13 + 3; %высота относительно уровня моря
    psi0 = atan(17.5/16.5);
    phi0 = (11+32/60+13.47/3600)/180*pi;
    lam0 = (104+50/60+03.58/3600)/180*pi;
    rn = 6378137*(1 - e2)/((1 - e2*(sin(phi0))^2)^(3/2));
    re = 6378137/((1 - e2*(sin(phi0))^2)^(1/2));
    
    T = [0 : DT : Lp/Vot]';
    gamma = T*0;
    theta = T*0;
    psi = psi0 + 0*T;
    H = H0 + T*0; %высота в каждый момент времени матрица размера как и Т ,но с числами H0
    Hsh(1, 1) = H0;
    x = T/Lp*Vot;
    x_ = 1/(T(end, 1) - T(1, 1));
    Ver = ((V0 - Vot)/2*cos(pi*x) + (V0 + Vot)/2)*sin(psi0);
    Vnr = ((V0 - Vot)/2*cos(pi*x) + (V0 + Vot)/2)*cos(psi0);
    V = [Ver, Vnr, T*0];
    Vdifer = -x_*pi*(V0 - Vot)/2*sin(pi*x)*sin(psi0);
    Vdifnr = -x_*pi*(V0 - Vot)/2*sin(pi*x)*cos(psi0);
    Vdif = [Vdifer, Vdifnr, T*0];
    Sgamma = [1, 0, 0; 0, 0, 1; 0, -1, 0];
    Spsi = [sin(psi0), cos(psi0), 0; -cos(psi0), sin(psi0), 0; 0, 0, 1];
    Sr = Sgamma*Spsi;
    phi(1, 1) = phi0;
    phish(1, 1) = phi0;
    lam(1, 1) = lam0;
    lamsh(1, 1) = lam0;
    G(1, :) = [0, 0, -(Ge*(cos(phi0))^2 + Gp*(1 - e2)^(1/2)*(sin(phi0))^2)/((1 - e2*(sin(phi0))^2)^(1/2))];
    Fx(:, 1) = -(G(1, :)');
    Fz(:, 1) = Sr*(-(G(1, :)'));
    Vsh = [T*0 T*0 T*0];
    Gsh(1, :) = [0, 0, -(Ge*(cos(phi0))^2 + Gp*(1 - e2)^(1/2)*(sin(phi0))^2)/(1 - e2*(sin(phi0))^2)^(1/2)];
    for i = 2 : size(T, 1)
        phi(i, 1) = phi(i - 1, 1) + Vnr(i, 1)*DT/(rn + H0);
        lam(i, 1) = lam(i - 1, 1) + Ver(i, 1)*DT/((re + H0)*cos(phi(i, 1)));
        G(i, :) = [0, 0, -(Ge*(cos(phi(i, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phi(i, 1)))^2)/((1 - e2*(sin(phi(i, 1)))^2)^(1/2))];
        Omega = [-V(i, 2)/(rn + H0), V(i, 1)/(re + H0), V(i, 1)*tan(phi(i, 1))/(re + H0)];
        Ux = [0, u*cos(phi(i, 1)), u*sin(phi(i, 1))];
        Omega_ = [0, Omega(3), -Omega(2); -Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
        Ux_ = [0, Ux(3), -Ux(2); -Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
        Fx(:, i) = (Vdif(i, :)' - (Omega_ + 2*Ux_)*V(i, :)' - G(i, :)');
        Fz(:, i) = Sr*(Vdif(i, :)' - (Omega_ + 2*Ux_)*V(i, :)' - G(i, :)');
        
        phish(i, 1) = phish(i - 1, 1) + Vsh(i - 1, 2)/(rn + Hsh(i - 1, 1))*DT;
        lamsh(i, 1) = lamsh(i - 1, 1) + Vsh(i - 1, 1)/(re + Hsh(i - 1, 1))/cos(phish(i - 1, 1))*DT;
        Hsh(i, 1) = Hsh(i - 1, 1) + Vsh(i - 1, 3)*DT;
        Gsh(i, :) = [0, 0, -(Ge*(cos(phish(i - 1, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phish(i - 1, 1)))^2)/(1 - e2*(sin(phish(i - 1, 1)))^2)^(1/2)];
        Omegash = [-Vsh(i - 1, 2)/(rn + Hsh(i - 1, 1)), Vsh(i - 1, 1)/(re + Hsh(i - 1, 1)), Vsh(i - 1, 1)*tan(phish(i - 1, 1))/(re + Hsh(i - 1, 1))];
        Omega_sh = [0, Omegash(3), -Omegash(2); -Omegash(3), 0, Omegash(1); Omegash(2), -Omegash(1), 0];
        Uxsh = [0, u*cos(phish(i - 1, 1)), u*sin(phish(i - 1, 1))];
        Ux_sh = [0, Uxsh(3), -Uxsh(2); -Uxsh(3), 0, Uxsh(1); Uxsh(2), -Uxsh(1), 0];
        o = (((Omega_sh + 2*Ux_sh)*Vsh(i - 1, :)')' + Gsh(i - 1, :) + (Sr'*Fz(:, i - 1))');
        Vsh(i, :) = Vsh(i - 1, :) + (((Omega_sh + 2*Ux_sh)*Vsh(i - 1, :)')' + Gsh(i - 1, :) + (Sr'*Fz(:, i - 1))')*DT;
    end
    FZr = [T, Fz'];
    flight0 = [T lam*180/pi phi*180/pi H gamma theta psi];
    
    insr = [T lamsh*180/pi phish*180/pi Hsh];
%  plot(insr(:, 1), Vsh(:, 3)); grid on;
    
    T1 = [T(end,1) : DT : T(end,1) + 15]';  
    x = (T1 - T(end,1))/(T1(1,1) - T1(end, 1));
    x_ = 1/(T1(1,1) - T1(end, 1));
    H1 = H0 + Hb;
    gammaot = T1*0;
    psiot = psi0 + 0*T1;
    v0 = 0;
    v1 = (H1 - H0)/5;
    theta1 = 15/180*pi;
    phi1 = phi(end, 1);
    lam1 = lam(end, 1);
    Veot = ((Vot - Vb)/2*cos(pi*x) + (Vb + Vot)/2)*sin(psi0);
    Vnot = ((Vot - Vb)/2*cos(pi*x) + (Vb + Vot)/2)*cos(psi0);
    Vupot = (v0 - v1)/2*cos(pi*x.^5) + (v0 + v1)/2;
    V2 = [Veot, Vnot, Vupot];
    V2difer = -x_*pi*(Vot - Vb)/2*sin(pi*x)*sin(psi0);
    V2difnr = -x_*pi*(Vot - Vb)/2*sin(pi*x)*cos(psi0);
    V2difup = -x_*pi*(5*x.^4)*(v0 - v1)/2.*sin(pi*(x).^(5));
    V2dif = [V2difer, V2difnr, V2difup];
    
    thetaot = - theta1/2*cos(pi*x) + theta1/2;
    phiot(1, 1) = phi1;
    lamot(1, 1) = lam1;
    Hot(1, 1) = H0;
    Got(1, :) = [0, 0, -(Ge*(cos(phiot(1, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiot(1, 1)))^2)/(1 - e2*(sin(phiot(1, 1)))^2)^(1/2)];
    Omega = [-V2(1, 2)/(rn + H0), V2(1, 1)/(re + H0), V2(1, 1)*tan(phiot(1, 1))/(re + Hot(1, 1))];
    Ux = [0, u*cos(phiot(1, 1)), u*sin(phiot(1, 1))];
    Omega_ = [0, Omega(3), -Omega(2); -Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
    Ux_ = [0, Ux(3), -Ux(2); -Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
    Fzot(:, 1) = Sr*(V2dif(1, :)' - (Omega_ + 2*Ux_)*V2(1, :)' - Got(1, :)');
    
    Votsh = [T1*0 T1*0 T1*0];
    phiotsh(1, 1) = phish(end, 1) + Vsh(end, 2)/(rn + Hsh(end, 1))*DT;
    lamotsh(1, 1) = lamsh(end, 1) + Vsh(end, 1)/(re + Hsh(end, 1))/cos(phish(end, 1))*DT;
    Hotsh(1, 1) = Hsh(end, 1) + Vsh(end, 3)*DT;
    Gotsh(1, :) = [0, 0, -(Ge*(cos(phish(end, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phish(end, 1)))^2)/(1 - e2*(sin(phish(end, 1)))^2)^(1/2)];
    Votsh(1, :) = Vsh(end, :) + DT*((Omega_sh + 2*Ux_sh)*Vsh(end, :)')' + Gsh(end, :) + (Sr'*Fz(:, end))';
    for i = 2 : size(T1, 1)
        Hot(i, 1) = Hot(i - 1, 1) + Vupot(i, 1)*DT;
        phiot(i, 1) = phiot(i - 1, 1) + Vnot(i, 1)*DT/(rn + Hot(i, 1));
        lamot(i, 1) = lamot(i - 1, 1) + Veot(i, 1)*DT/((re + Hot(i, 1))*cos(phiot(i, 1)));
        Got(i, :) = [0, 0, -(Ge*(cos(phiot(i, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiot(i, 1)))^2)/(1 - e2*(sin(phiot(i, 1)))^2)^(1/2)];
        Omega = [-V2(i, 2)/(rn + Hot(i, 1)), V2(i, 1)/(re + Hot(i, 1)), V2(i, 1)*tan(phiot(i, 1))/(re + Hot(i, 1))];
        Ux = [0, u*cos(phiot(i, 1)), u*sin(phiot(i, 1))];
        Omega_ = [0, Omega(3), -Omega(2); -Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
        Ux_ = [0, Ux(3), -Ux(2); -Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
        Stheta = [cos(thetaot(i, 1)), 0, sin(thetaot(i, 1)); 0, 1, 0; -sin(thetaot(i, 1)), 0, cos(thetaot(i, 1))];
        S = Sgamma*Stheta*Spsi;
        Fzot(:, i) = S*(V2dif(i, :)' - (Omega_ + 2*Ux_)*V2(i, :)' - Got(i, :)');
        
        phiotsh(i, 1) = phiotsh(i - 1, 1) + Votsh(i - 1, 2)/(rn + Hotsh(i - 1, 1))*DT;
        lamotsh(i, 1) = lamotsh(i - 1, 1) + Votsh(i - 1, 1)/(re + Hotsh(i - 1, 1))/cos(phiotsh(i - 1, 1))*DT;
        Hotsh(i, 1) = Hotsh(i - 1, 1) + Votsh(i - 1, 3)*DT;
        Gotsh(i, :) = [0, 0, -(Ge*(cos(phiotsh(i - 1, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiotsh(i - 1, 1)))^2)/(1 - e2*(sin(phiotsh(i - 1, 1)))^2)^(1/2)];
        Omegash = [-Votsh(i - 1, 2)/(rn + Hotsh(i - 1, 1)), Votsh(i - 1, 1)/(re + Hotsh(i - 1, 1)), Votsh(i - 1, 1)*tan(phiotsh(i - 1, 1))/(re + Hotsh(i - 1, 1))];
        Omega_sh = [0, Omegash(3), -Omegash(2); -Omegash(3), 0, Omegash(1); Omegash(2), -Omegash(1), 0];
        Uxsh = [0, u*cos(phiotsh(i - 1, 1)), u*sin(phiotsh(i - 1, 1))];
        Ux_sh = [0, Uxsh(3), -Uxsh(2); -Uxsh(3), 0, Uxsh(1); Uxsh(2), -Uxsh(1), 0];
        Votsh(i, :) = Votsh(i - 1, :) + (((Omega_sh + 2*Ux_sh)*Votsh(i - 1, :)')' + Gotsh(i - 1, :) + (S'*Fzot(:, i - 1))')*DT;
    end
    FZot = [T1, Fzot'];
    flight1 = [T1 lamot*180/pi phiot*180/pi Hot gammaot thetaot psiot];
    insot = [T1 lamotsh*180/pi phiotsh*180/pi Hotsh];
    
    T2 = [ T1(end,1): DT : T1(end,1) + 700]';
    Vkr = 296/3600; 
    x = (T2 - T2(1, 1))/(T2(end, 1) - T2(1, 1)); 
    x_ = 1/(T2(end, 1) - T2(1, 1));
    Veup = ((Vb - Vkr)/2*cos(pi*x) + (Vkr + Vb)/2)*sin(psi0); 
    Vnup = ((Vb - Vkr)/2*cos(pi*x) + (Vkr + Vb)/2)*cos(psi0); 
    Vup = (Vupot(end, 1) - Vkr)/2*(cos(pi*(1 - (1 - x).^2))) + (Vupot(end, 1) + Vkr)/2; 
    V3 = [Veup, Vnup, Vup];
    V3difer = -x_*pi*(Vb - Vkr)/2*sin(pi*x)*sin(psi0);
    V3difnr = -x_*pi*(Vb - Vkr)/2*sin(pi*x)*cos(psi0);
    V3difup = -x_*pi*(1 - x)*(Vupot(end, 1) - Vkr)/2.*sin(pi*x);
    V3dif = [V3difer, V3difnr, V3difup];
   
    gammaup = 0*T2;
    thetaup = thetaot(end, 1)/2*cos(pi*x) + thetaot(end, 1)/2; 
    psiup = psi0 + 0*T2; 
    phiup(1, 1) = phiot(end, 1); 
    lamup(1, 1) = lamot(end, 1); 
    Hup(1, 1) = Hot(end, 1);
 
    Gup(1, :) = [0, 0, -(Ge*(cos(phiot(1, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiot(1, 1)))^2)/(1 - e2*(sin(phiot(1, 1)))^2)^(1/2)];
    Omega = [-V3(1, 2)/(rn + Hup(1, 1)), V3(1, 1)/(re + Hup(1, 1)), V3(1, 1)*tan(phiot(1, 1))/(re + Hup(1, 1))];
    Ux = [0, u*cos(phiup(1, 1)), u*sin(phiup(1, 1))];
    Omega_ = [0, Omega(3), -Omega(2); -Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
    Ux_ = [0, Ux(3), -Ux(2); -Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
    Stheta = [cos(thetaup(1, 1)), 0, sin(thetaup(1, 1)); 0, 1, 0; -sin(thetaup(1, 1)), 0, cos(thetaup(1, 1))];
    S = Sgamma*Stheta*Spsi;
    Fzup(:, 1) = S*(V3dif(1, :)' - (Omega_ + 2*Ux_)*V3(1, :)' - Gup(1, :)');
    
    Vupsh = [T2*0 T2*0 T2*0];
    phiupsh(1, 1) = phiotsh(end, 1) + Votsh(end, 2)/(rn + Hotsh(end, 1))*DT;
    lamupsh(1, 1) = lamotsh(end, 1) + Votsh(end, 1)/(re + Hotsh(end, 1))/cos(phiotsh(end, 1))*DT;
    Hupsh(1, 1) = Hotsh(end, 1) + Votsh(end, 3)*DT;
    Gupsh(1, :) = [0, 0, -(Ge*(cos(phiotsh(end, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiotsh(end, 1)))^2)/(1 - e2*(sin(phiotsh(end, 1)))^2)^(1/2)];
    Vupsh(1, :) = Votsh(end, :) + DT*((Omega_sh + 2*Ux_sh)*Votsh(end, :)')' + Gotsh(end, :) + (S'*Fzot(:, end))';
    for i = 2 : size(T2, 1) 
         Hup(i, 1) = Hup(i - 1, 1) + Vup(i, 1)*DT; 
         phiup(i, 1) = phiup(i - 1, 1) + Vnup(i, 1)*DT/(rn + Hup(i, 1)); 
         lamup(i, 1) = lamup(i - 1, 1) + Veup(i, 1)*DT/((re + Hup(i, 1))*cos(phiup(i, 1)));
         Gup(i, :) = [0, 0, -(Ge*(cos(phiup(i, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiup(i, 1)))^2)/(1 - e2*(sin(phiup(i, 1)))^2)^(1/2)];
         Omega = [-V3(i, 2)/(rn + Hup(i, 1)), V3(i, 1)/(re + Hup(i, 1)), V3(i, 1)*tan(phiup(i, 1))/(re + Hup(i, 1))];
         Ux = [0, u*cos(phiup(i, 1)), u*sin(phiup(i, 1))];
         Omega_ = [0, Omega(3), -Omega(2); -Omega(3), 0, Omega(1); Omega(2), -Omega(1), 0];
         Ux_ = [0, Ux(3), -Ux(2); -Ux(3), 0, Ux(1); Ux(2), -Ux(1), 0];
         Stheta = [cos(thetaup(i, 1)), 0, sin(thetaup(i, 1)); 0, 1, 0; -sin(thetaup(i, 1)), 0, cos(thetaup(i, 1))];
         S = Sgamma*Stheta*Spsi;
         Fzup(:, i) = S*(V3dif(i, :)' - (Omega_ + 2*Ux_)*V3(i, :)' - Gup(i, :)');
         
        phiupsh(i, 1) = phiupsh(i - 1, 1) + Vupsh(i - 1, 2)/(rn + Hupsh(i - 1, 1))*DT;
        lamupsh(i, 1) = lamupsh(i - 1, 1) + Vupsh(i - 1, 1)/(re + Hupsh(i - 1, 1))/cos(phiupsh(i - 1, 1))*DT;
        Hupsh(i, 1) = Hupsh(i - 1, 1) + Vupsh(i - 1, 3)*DT;
        Gupsh(i, :) = [0, 0, -(Ge*(cos(phiupsh(i - 1, 1)))^2 + Gp*(1 - e2)^(1/2)*(sin(phiupsh(i - 1, 1)))^2)/(1 - e2*(sin(phiupsh(i - 1, 1)))^2)^(1/2)];
        Omegash = [-Vupsh(i - 1, 2)/(rn + Hupsh(i - 1, 1)), Vupsh(i - 1, 1)/(re + Hupsh(i - 1, 1)), Vupsh(i - 1, 1)*tan(phiupsh(i - 1, 1))/(re + Hupsh(i - 1, 1))];
        Omega_sh = [0, Omegash(3), -Omegash(2); -Omegash(3), 0, Omegash(1); Omegash(2), -Omegash(1), 0];
        Uxsh = [0, u*cos(phiupsh(i - 1, 1)), u*sin(phiupsh(i - 1, 1))];
        Ux_sh = [0, Uxsh(3), -Uxsh(2); -Uxsh(3), 0, Uxsh(1); Uxsh(2), -Uxsh(1), 0];
        o(i, :) = (((Omega_sh + 2*Ux_sh)*Vupsh(i - 1, :)')' + Gupsh(i - 1, :) + (S'*Fzup(:, i - 1))');
        Vupsh(i, :) = Vupsh(i - 1, :) + (((Omega_sh + 2*Ux_sh)*Vupsh(i - 1, :)')' + Gupsh(i - 1, :) + (S'*Fzup(:, i - 1))')*DT;
    end 
    FZup = [T2, Fzup'];
    FZ = [FZr; FZot; FZup];
    insup = [T2 lamupsh*180/pi phiupsh*180/pi Hupsh];

    ins = [insr; insot; insup];
    
    %save ins.txt ins �ascii;
   % tgeo2kml(ins(1:10:end,:));
    save fz.txt FZ -ascii;
    flight2 = [T2 lamup*180/pi phiup*180/pi Hup gammaup thetaup psiup]; 
    flight = [flight0 ; flight1; flight2 ];
     %plot(ins(:, 1), ins(:, 4)); grid on;
    save flight.txt flight -ascii;
    tgeo2kml(flight)
    
              