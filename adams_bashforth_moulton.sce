clc;
function  res = rungekutta4O(xn, yn, h, f )
    
    k1 = h*f(xn, yn)
    k2 = h*f(xn+h/2 , yn+k1/2)
    k3 = h*f(xn+h/2 , yn+k2/2)
    k4 = h*f(xn+h , yn+k3)
    
    res = yn + (k1 + 2*k2 + 2*k3 + k4)/6
endfunction

function res=adams(x0, y0, xn, precision, f)
    h=(xn-x0)/4;
    x1=x0+h;    x2=x1+h;    x3=x2+h;    x4=x3+h;
    y1=rungekutta4O(x0, y0, h, f)
    y2=rungekutta4O(x1, y1, h, f)
    y3=rungekutta4O(x2, y2, h, f)
    y0_=f(x0,y0);
    y1_=f(x1, y1);
    y2_=f(x2, y2);
    y3_=f(x3, y3);
    p4=y3+(55*y3_ - 59*y2_ + 37*y1_ - 9*y0_)*h/24;
    y4_=f(x4, p4);  y4=0;
    
    dif=200;
    while dif> (0.5 * 10^(-1*precision))
        p4=y4;
        y4=y3+(9*y4_ + 19*y3_ - 5*y2_ + y1_)*h/24;
        y4_=f(x4, y4);
        dif=abs(y4-p4)
    end
    
    res=y4
endfunction


deff('g=f(x,y)', 'g=y-x^2');
adams(0, 1, 0.4, 5, f)
