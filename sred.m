function osred

  global ep;
  
  e=2.7;
  ep = 0.01;
  dt = 1e-2;
  t = [0: dt: 20];
  y10 = 0;
  y20 = 0;
  y1 = [y10 y20];
  y2 = [y10 y20];
  for i = 2:size(t, 2)
    
    k1 = f(t(i - 1), y1(i - 1, :));
    k2 = f(t(i - 1) + dt/2, y1(i - 1, :) + dt/2*k1);
    k3 = f(t(i - 1) + dt/2, y1(i - 1, :) + dt/2*k2);
    k4 = f(t(i - 1) + dt, y1(i - 1, :) + dt*k3);
    y1(i, :) = y1(i - 1, :) + dt/6*(k1 + 2*k2 + 2*k3 + k4);
   
    
    y2(i, 1) = (e^(-0.005*t(i)))*(((-0.0001*e^(-0.005*t(i))*sin(2*t(i))+0.04*cos(2*t(i)))/16.0001)+0.98638-((-0.0001*e^(-0.005*t(i))*sin(2*t(i))+0.04*cos(2*t(i)))/16.0001)-0.5*e^(-0.005*t(i)));
    y2(i, 2) =-((e^(-0.005*t(i)))*(((-0.0001*e^(-0.005*t(i))*sin(2*t(i))+0.04*cos(2*t(i)))/16.0001)+0.98638-((-0.0001*e^(-0.005*t(i))*sin(2*t(i))+0.04*cos(2*t(i)))/16.0001)-0.5*e^(-0.005*t(i)))+sin(4*t(i))+0.5*sin(2*t(i)));
    
  end
  plot(t, y1(:, 1)-y2(:, 1));
  grid on;
  pause(2);
 % close all;
  plot(t, y1(:, 2), t, y2(:, 2),'r');
  grid on;
  
 end
 
 function dydt = f(t, y)
 
  global ep;
  dydt = [0.5*y(1)+y(2)+0.5*(sin(4*t)+0.5*sin(2*t)), y(1)+sin(4*t)+0.5*sin(2*t)+ep];
 
 end