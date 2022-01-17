function puan

  global mu;
  
  mu = 0.1;
  dt = 1e-3;
  t = [0: dt: 1];
  y10 = 0;
  y20 = 1;
  y1 = [y10 y20];
  y2 = [y10 y20];
  for i = 2:size(t, 2)
    
    k1 = f(t(i - 1), y1(i - 1, :));
    k2 = f(t(i - 1) + dt/2, y1(i - 1, :) + dt/2*k1);
    k3 = f(t(i - 1) + dt/2, y1(i - 1, :) + dt/2*k2);
    k4 = f(t(i - 1) + dt, y1(i - 1, :) + dt*k3);
    y1(i, :) = y1(i - 1, :) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    
    y2(i, 1) = t(i)^4/4;
    y2(i, 2) = -(t^9/144)+1;
  end
  
  plot(t, y1(:, 1), t, y2(:, 1));
  pause(2);
  close all;
  plot(t, y1(:, 2), t, y2(:, 2));;
  
 end
 
 function dydt = f(t, y)
 
  global mu;
  dydt = [t^3*(cos(mu*t)), -y(1)^2 + cos(mu*y(1))];
 
 end