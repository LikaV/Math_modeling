function psyche

  [x, y, z] = sphere(24);
  surf(120e3*x, 90e3*y, 70e3*z);
  daspect([1 1 1]);
   hold on;

  X_ = 500e3;
  T_ = 5e3;
  D=0;
 %s_all=[0.1 0.01 0.001]
  s_all=10.^(-[1:7]);
  s0 =1e-2/0.6;
  for k=1:7
  x0 = [-1/2 0 1/5];
  v0 = [1 0 0];
  y1 = [x0 v0];
  y2 = [x0 v0];
  s0 =1e-2/0.6;
  s=s_all(k);
  dt  = 1e-2;
  t = [0:dt:1]';

  
  for i  = 2:size(t,1)
  	
   
    i = i + 1;
		t(i, 1) = t(i - 1, 1) + dt;
	
    k1 = f1(t(i - 1), y1,s);
    k2 = f1(t(i - 1) + dt/2, y1 + dt/2*k1,s);
    k3 = f1(t(i - 1) + dt/2, y1 + dt/2*k2,s);
    k4 = f1(t(i - 1) + dt, y1 + dt*k3,s);
    y1 = y1 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    
    k1 = f2(t(i - 1), y2);
    k2 = f2(t(i - 1) + dt/2, y2 + dt/2*k1);
    k3 = f2(t(i - 1) + dt/2, y2 + dt/2*k2);
    k4 = f2(t(i - 1) + dt, y2 + dt*k3);
    y2 = y2 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    
    d=(y1-y2)*(y1-y2)';
    if d > D
        D=d;
    end
    D_all(k)=D;
    
    plot3(X_*y1(1), X_*y1(2), X_*y1(3), 'r.');
    plot3(X_*y2(1), X_*y2(2), X_*y2(3), 'b.');
    pause(dt);
  end
     hold off;
     
     loglog(s_all, D_all,'*-');
     xlabel('s');
     ylabel('\Delta');
     grid on;
  end
  
end
 function dydt = f1(t, y, s)
 
  v =  y(4:6);
  dydt(1:3) = v;
  x = y(1:3);
  dydt(4:6) = -0.3*x/norm(x)^3 + s*cross(v, b(x));
  
 end

  function dydt = f2(t, y)
 
  v =  y(4:6);
  dydt(1:3) = v;
  x = y(1:3);
  dydt(4:6) = -0.3*x/norm(x)^3 + 0*cross(v, b(x));
  
 end
 
 function bv = b(x) 
  
  phi = atan(x(3)/norm(x(1:2)));
  r = norm(x);
  r0 = 180/500; %R0 nosize
  
  y = [-x(1)*sin(phi)/norm(x(1:2)), -x(2)*sin(phi)/norm(x(1:2)), cos(phi)];
  z = x/r;
  
  bN = -(r0/r)^3*(-cos(phi) - 3/2*(1/10)*(r0/r)*sin(2*phi)); % единицу измерения гамма берем как модуль
  bUP = (r0/r)^3*(-2*sin(phi) - 3/2*(r0/r)*(3*sin(phi)^2 - 1));
  
  bv = bN*y + bUP*z;
  
 end
  