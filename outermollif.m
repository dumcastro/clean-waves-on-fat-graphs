function P_ep = outermollif(P,ep,angles,widths)

    %% Process polygon vertices
    ver = vertex(P);
    
    
    % Apply small extension to avoid singularities
    if size(P) == [9 1]
    %
    thet = [0, angles(2), angles(2)+3*pi/2, angles(2), angles(3),...
        angles(3) + pi/2, angles(3)+pi/2, angles(3)];
    thet = pi + thet;

    ver(1) = ver(1) -ep*1i;
    ver(9) = ver(9) + ep*1i;

    ver(3) = ver(3) + ep*exp(1i*thet(3));
    ver(7) = ver(7) + ep*exp(1i*thet(7));

    lmbd = (-ep-imag(ver(3)))/sin(thet(2));

    ver(2) = ver(3) + lmbd*exp(1i*thet(2));

    lmbd = (widths(1)+ep-imag(ver(7)))/sin(thet(8));

    ver(8) = ver(7) + lmbd*exp(1i*thet(8));


    %{
    for j = 2:8
        
        A = [cos(thet(j-1)), - cos(thet(j-1)+(thet(j)-thet(j-1))/2) ;
            sin(thet(j-1)), - sin(thet(j-1)+(thet(j)-thet(j-1))/2)];

        b = [real(ver(j)-ver(j-1));imag(ver(j)-ver(j-1))];

        lmdb = A\b;

        ver(j) = ver(j-1) + lmdb(1)*exp(1i*thet(j-1));

    end
    %}

    ver = ver + ep*1i;

    %}
    else
    ver(2) = ver(2) - ep*1i;
    ver(5) = ver(5) + ep*1i;
    end

    P_ep = polygon(ver)

end