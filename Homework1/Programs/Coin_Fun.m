function eta  = Coin_Fun( xi )
% this function tranfer one stochastic varibale xi to another stochastic varibale eta  
% Which means the sturcture of two side of coin
if xi < 0.5
    eta = 0;
else
    eta = 1;
end

end

