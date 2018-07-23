function [out] = translate (pos, delta)
D = ones(size(pos),'single').*delta;
out = pos + D;
end