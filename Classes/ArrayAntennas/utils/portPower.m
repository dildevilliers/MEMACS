function P = portPower(x)
% Returns the power at the antenna ports
% x is the getPortSignal matrix [Nant x Nsamp]
% Of course not very sensible for 1 or 2 bit sampling...
P = sum(abs(x).^2,2)./size(x,2);
end