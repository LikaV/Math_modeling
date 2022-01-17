function pendulum

  global G L

	PHI0 = pi/3; %rad
	G = 1.6; %m/s2
	L = 0.1; % m
	OMEGA0 = 0; % rad/s

	T1 = 50; % sqrt(8*PHI0*L/G/sin(PHI0));
	DT = 1e-2; % 1e-3*T1;

	Y = [PHI0 OMEGA0];
	T = 0;
	i = 1; % N �����
	initdraw([0 0]', Y', L, L/10);
	while T(i,1) < T1
	
		i = i + 1;
		T(i,1) = T(i-1,1) + DT;
	
		%PHI = Y(i-1,1);
		%OMEGA = Y(i-1,2);
		%Y(i,:) = Y(i-1,:) + [OMEGA, -G/L*sin(PHI)]*DT;
    k1 = F(T(i - 1, 1), Y(i -1, :));
    k2 = F(T(i - 1, 1) + DT/2, Y(i -1, :) + DT/2*k1);
    k3 = F(T(i - 1, 1) + DT/2, Y(i -1, :) + DT/2*k2);
    k4 = F(T(i - 1, 1) + DT, Y(i - 1, :) + DT*k3);
    Y(i,:) = Y(i - 1, :) + DT/6*(k1 + 2*k2 + 2*k3 + k4);
	
		if rem(i,10) == 0
			stepdraw([0 0]', Y(i,:)', L, L/10);
			pause(DT);
		end
	
	end

	plot(T,Y(:,1));
	grid on;
	pause(1);
	
end

function dYdT = F(T, Y)
  
  global G L
  
  dYdT = [Y(2), -G/L*sin(Y(1))]; 
  
end


function initdraw(H, X, L, R)
    
    global hl;    global hb
    global hv
    global hh
    
    f = X(1,1);
    v = X(2,1);
    
    close all;
    hl = line([H(1,1) H(1,1)+L*sin(f)],[H(2,1) H(2,1)-L*cos(f)],'LineWidth',2);
    %hold on;
    hb = rectangle('Position',[H(1,1)+L*sin(f)-R, H(2,1)-L*cos(f)-R, 2*R, 2*R], ...
            'Curvature', 1, 'FaceColor', 'g');
    hv = line(  [H(1,1)+(L+R)*sin(f) H(1,1)+(L+R)*sin(f)+v*cos(f)], ...
                [H(2,1)-(L+R)*cos(f) H(2,1)-(L+R)*cos(f)+v*sin(f)], ...
                'Color','m');
    hh = rectangle('Position', [H(1,1)-L/50, H(2,1)-L/50, L/25, L/25],...
        'Curvature', 1, 'FaceColor','k');
    %hold off;
    
    axis([-L L -L L]*3);
    daspect([1 1 1]);
    grid on;
    
end

function stepdraw(H, X, L, R)

    global hl
    global hb
    global hv
    global hh
    
    f = X(1,1);
    v = X(2,1);
    
    set(hl, 'XData', [H(1,1) H(1,1)+L*sin(X(1,1))], ...
            'YData', [H(2,1) H(2,1)-L*cos(X(1,1))]);
    set(hv, 'XData', [H(1,1)+(L+R)*sin(f) H(1,1)+(L+R)*sin(f)+v*cos(f)], ...
            'YData', [H(2,1)-(L+R)*cos(f) H(2,1)-(L+R)*cos(f)+v*sin(f)]);
    
    set(hb, 'Position',...
        [H(1,1)+L*sin(X(1,1))-R, H(2,1)-L*cos(X(1,1))-R, 2*R, 2*R]);
    
    set(hh, 'Position', [H(1,1)-L/50, H(2,1)-L/50, L/25, L/25]);
    
    drawnow;

end
