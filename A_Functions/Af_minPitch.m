function minPitch = Af_minPitch(u_bar,x0,x1,fine_pitch,a)

% a  = 10;
% x0 = 12;
% x1 = 20;

if u_bar < x0
    minPitch = fine_pitch;
elseif u_bar > x1
    minPitch = a + fine_pitch;
else
    a3 = 2/(x0-x1)^3;
    a2 = -3*(x0+x1)/(x0-x1)^3;
    a1 = 6*x1*x0/(x0-x1)^3;
    a0 = (x0-3*x1)*x0^2/(x0-x1)^3;
    
    sigma = [a3,a2,a1,a0] * [u_bar.^3;u_bar.^2;u_bar;ones(size(u_bar))];
    
    minPitch = a*sigma + fine_pitch;
end