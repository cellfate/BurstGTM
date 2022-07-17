function [theta,theta0] = burstInference(vals)

theta0 = momentInference(vals);
if isempty(theta0) || any(theta0 <= 0)
    theta0 = [10,10,20];
end

theta = MLInference(vals,theta0);

end
