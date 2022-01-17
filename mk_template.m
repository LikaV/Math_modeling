function mk

    loadimg;

    t = 1; % turn number

    H1 = 100; H2 = 100; % player's health
	dH1 = 0; dH2 = 0;   % player's health decrease

    while 1 ...
		
		showfight(1,1, H1, H2, dH1, dH2);
        pause();
		
		a1 = 1; % player 1 action
		a2 = 1; % player 2 action
        
		...
		
        showfight(a1,a2,H1,H2,dH1,dH2);
        pause();
    
        t = t + 1;
    
	end

end

function showfight(a1,a2,H1,H2,dH1,dH2)

    global IMG;
    
    if H1 > 0
        i1 = IMG{1,a1};
    else
        i1 = IMG{1,end};
        disp('Sub-zero wins');
    end
    if H2 > 0
        i2 = IMG{2,a2};
    else
        i2 = IMG{2,end};
        disp('Scorpion wins');
    end
    n1 = size(i1,1);
    m1 = size(i1,2);
    n2 = size(i2,1);
    m2 = size(i2,2);
    
	for k = 1:3
		i(:,:,k) = [ [i1(:,:,k); 255*ones(n2-n1,m1)], ...
			[i2(:,:,k); 255*ones(n1-n2,m2)]];
	end
    n = max([n1 n2]);
    m = m1 + m2;

    imshow(i, 'XData', [-m/2 +m/2]*10, 'YData', [-n/2 +n/2]*10);
    axis image
    truesize
    text(-m/2-1000, -20, num2str(max(H1,0)));
    text(+m/2+1000, -20, num2str(max(H2,0)));
	text(-m/2-1000, +120, num2str(dH1));
    text(+m/2+1000, +120, num2str(dH2));

end

function loadimg

    global IMG
    
    IMG = cell(2,1);
    
	for i = 1:6

		IMG(1,i) = {imread(['picmk\a' num2str(i) '.jpg'])};
		IMG(2,i) = {imread(['picmk\b' num2str(i) '.jpg'])};

	end

end