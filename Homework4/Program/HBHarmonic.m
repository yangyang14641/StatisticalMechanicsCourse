function dy = HBHarmonic(t,y)
% Hamilton equations for heat bath harmonic oscillator
global Ms a b

dy = zeros(3,1);    % a column vector

dy(1) = -a*y(2) - y(1)*y(3);
dy(2) = b * y(1);
dy(3) = 1/Ms * (y(1)^2 - 1);
end

