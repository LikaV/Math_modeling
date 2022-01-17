function gravity

    global Xe Xv Xm
    global Re Rv Rm Lm 
    global Me Mv Mm G
	
  Me = 6e24; %кг
	Re = 6e6; 
  Mv = 5e24;
  Rv = 6e6;
  Mm = 7e23; 
  Rm = 2e6; 
  Lm = 4e8;
  Um = 1e3;
  G = 20/3*1e-11;
  T1 = 6e5;
  DT = 1e-3*T1;
  T = [0: DT: 5*T1]';
  N = size(T, 1);
  Y = [0 0 0  0 0 0  0 Lm 0  0.9e3 0 0  Lm 0 0  0 Um 0];
  
	Xe = Y(1:3);
  Xv = Y(7:9);
  Xm = Y(13:15);
  init_animation;   
    for i = 2:N  
     
      k1 = F(T(i - 1, 1), Y(i -1, :));
      k2 = F(T(i - 1, 1) + DT/2, Y(i -1, :) + DT/2*k1);
      k3 = F(T(i - 1, 1) + DT/2, Y(i -1, :) + DT/2*k2);
      k4 = F(T(i - 1, 1) + DT, Y(i - 1, :) + DT*k3);
      Y(i,:) = Y(i - 1, :) + DT/6*(k1 + 2*k2 + 2*k3 + k4);
        
     Xm = Y(i, 13:15);
     Xv = Y(i, 7:9);
     Xe = Y(i, 1:3);
       if rem(i,10) == 0
  	  	  step_animation; 
		  	  pause(DT/T1); 
	     end
        
       if norm(Xe - Xm) < Re + Rm
          disp('BOOOM')
          break
       end
       if norm(Xv - Xm) < Rv + Rm
          disp('BOOOM')
          Y(i, 10) = (Mm/(Mm + Mv)*(Y(i, 10))^2 + Mv/(Mm + Mv)*(Y(i, 16))^2)^(1/2);
          Y(i, 11) = (Mm/(Mm + Mv)*(Y(i, 11))^2 + Mv/(Mm + Mv)*(Y(i, 17))^2)^(1/2);
          disp(Y(i, :))
          Mv = Mv + Mm;
          Rm = 0; Mm = 0; Y(i, 16:18) = 0; Y(i, 12:15) = 0;
       end
       if norm(Xe - Xv) < Re + Rv
          disp('BOOOM')
          break
       end
        R (i, 1) = norm(Xm - Xv);
    end
    close all;
    plot(R);
    grid on
    
    
end

function dYdT = F(T, Y)

     global Me Mv Mm G
     dYdT(1:3)  = Y(4:6);
     dYdT(7:9)  = Y(10:12);
     dYdT(13:15)  = Y(16:18);
     
     Xm = Y(13:15);
     Xv = Y(7:9);
     Xe = Y(1:3);
     Rve =norm(Xv  - Xe); 
     Rme = norm(Xm  - Xe);
     Rmv = norm(Xm  - Xv);
     dYdT(4:6)  = G*Mm/Rme^2*(Xm - Xe)/Rme +  G*Mv/Rve^2*(Xv - Xe)/Rve;
     dYdT(10:12)  = G*Mm/Rmv^2*(Xm - Xv)/Rmv - G*Me/Rve^2*(Xv - Xe)/Rve;
     dYdT(16:18)  = - G*Me/Rme^2*(Xm - Xe)/Rme - G*Mv/Rmv^2*(Xm - Xv)/Rmv;

end     

function init_animation

    global Xe Xv Xm 
    global Re Rv Rm Lm
    global he hv hl
    global xs ys zs
    
    [xs,ys,zs] = sphere(12);
    
    close all;
    figure;
    hold off;
    colormap('winter');
    he = surf(Re*xs+Xe(1), Re*ys+Xe(2), Re*zs+Xe(3));
    hold on;
    hv = surf(Rv*xs+Xv(1), Rv*ys+Xv(2), Rv*zs+Xv(3));
    hl = surf(Rm*xs+Xm(1), Rm*ys+Xm(2), Rm*zs+Xm(3));
    
    grid on;
    daspect([1 1 1]);
    axis([-Lm 3*Lm -3*Lm Lm -Re Re]);
                    
    set(gcf, 'units','normalized', 'OuterPosition', [0 0 1 1]);

end

function step_animation

    global Xe Xv Xm
    global Re Rv Rm
    global he hv hl
    global xs ys zs
    
    set(he, ...
        'XData', Re*xs+Xe(1), ...
        'YData', Re*ys+Xe(2), ...
        'ZData', Re*zs+Xe(3));
    set(hv, ...
        'XData', Rv*xs+Xv(1), ...
        'YData', Rv*ys+Xv(2), ...
        'ZData', Rv*zs+Xv(3));
    set(hl, ...
        'XData', Rm*xs+Xm(1), ...
        'YData', Rm*ys+Xm(2), ...
        'ZData', Rm*zs+Xm(3));
    plot3(Xe(1), Xe(2), Xe(3), 'g.');
    plot3(Xv(1), Xv(2), Xv(3), 'b.');
    plot3(Xm(1), Xm(2), Xm(3), 'r.');

end
    