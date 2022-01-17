function two

m_all= 10.^(-[1:7])';
for k = 1:size(m_all,1)
x0 = [0];
v0 = [1];
y1 = [x0 v0];
y2 = [x0 v0];
dt = 1e-2;
t = [0:dt:1]';
m = m_all(k);
D = 0;
for i = 2:40
      k1 = f1( t(i-1), y1, m );
      k2 = f1( t(i-1) + dt/2, y1 + dt/2*k1, m);
      k3 = f1( t(i-1) + dt/2, y1 + dt/2*k2, m);
      k4 = f1( t(i-1) + dt, y1 + dt*k3, m);        
      y1 = y1 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
      

      k1 = f2( t(i-1), y2 );
      k2 = f2( t(i-1) + dt/2, y2 + dt/2*k1);
      k3 = f2( t(i-1) + dt/2, y2 + dt/2*k2);
      k4 = f2( t(i-1) + dt, y2 + dt*k3);    
      y2 = y2 + dt/6*(k1 + 2*k2 + 2*k3 + k4);
      
      d = (y1-y2)*(y1-y2)';
      if d>D
        D = d;
      end
      pause(dt/10)
end
D_all(k,1) = D;
end

close all;
loglog(m_all, D_all, '*-');
xlabel('\mu');
ylabel('\Delta');
grid on;
hold off;

end

function dydt = f1(t,y, m)

x = y(1);
v = y(2);
dydt(1) = v;
dydt(2) = (2 - exp(-t))/(-m*t);


end


function dydt = f2(t,y)

x = y(1);
v = y(2);
dydt(1) = v;
dydt(2) =  (2 - exp(-t))/(-0*t);


end