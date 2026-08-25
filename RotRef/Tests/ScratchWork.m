
test = @(x) soft_weight(x, 15, 2);

figure;
fplot(test,[10 20]);

function w = soft_weight(x, threshold, delta)
%Cubic Hermite smoothstep:
%   w = 0           for x <= threshold - delta
%   w = 1           for x >= threshold + delta
%   w = smooth      for threshold - delta < x < threshold + delta
%with w'=0 at both endpoints so the composition w(Dim(theta)) is C^1 in
%theta away from a zero-measure set.
if delta <= 0
    %Degenerate to hard indicator
    w = double(x > threshold);
    return;
end
t = (x - (threshold - delta)) / (2*delta);
t = max(0, min(1, t));
w = t .* t .* (3 - 2*t);
end