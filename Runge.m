function arhimed

global G P

M = 17; %кг
R = 1100; %кг/м^3
S = 1; %м
P = R*S/M; %1/м
G = 9.8; %м/с^2

DT = 0.1;
T = 0;
T1 = 5;
Z = 0;
V = 0;
Y = [Z, V];
i = 1;

  while T(i, 1) < T1
	
		i = i + 1;
		T(i, 1) = T(i - 1, 1) + DT;
	
    k1 = F(T(i - 1, 1), Y(i - 1, :));
    k2 = F(T(i - 1, 1) + DT/2, Y(i - 1, :) + DT/2*k1);
    k3 = F(T(i - 1, 1) + DT/2, Y(i - 1, :) + DT/2*k2);
    k4 = F(T(i - 1, 1) + DT, Y(i - 1, :) + DT*k3);
    Y(i, :) = Y(i - 1, :) + DT/6*(k1 + 2*k2 + 2*k3 + k4);
	  
	end
  k = Y(i - 1, 1)*P %Коэффициент в формуле зависимости
  plot(T(:, 1), Y(:, 1));
  grid on;
 

end
 
function dYdT = F(T, Y)
  
  global G P
  
  dYdT = [Y(2), -1*G*(1 + P*Y(1))]; 
  
end