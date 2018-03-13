function dy = Harmonic(t,y)
% Hamilton equations for harmonic oscillator
global a b

dy = zeros(2,1);    % a column vector

dy(1) = -a*y(2);
dy(2) = b * y(1);
end