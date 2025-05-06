% Function to create an interpolated model for a specific day
function model_fun = createInterpolatedModel(day, windowCenters, modelFunctions)
    % Find the two closest window centers
    [dists, idx] = sort(abs(windowCenters - day));
    
    % If we're exactly at a window center or only have one window, return that model
    if dists(1) == 0 || length(windowCenters) == 1
        model_fun = modelFunctions{idx(1)};
        return;
    end
    
    % Get the two closest models
    model1 = modelFunctions{idx(1)};
    model2 = modelFunctions{idx(2)};
    
    % Calculate interpolation weights
    day1 = windowCenters(idx(1));
    day2 = windowCenters(idx(2));
    
    % Ensure day is between day1 and day2
    if day < min(day1, day2) || day > max(day1, day2)
        % Day is outside the range, use the nearest model
        model_fun = modelFunctions{idx(1)};
        return;
    end
    
    % Linear interpolation weight
    if day2 ~= day1  % Avoid division by zero
        w = (day - day1) / (day2 - day1);
    else
        w = 0.5;
    end
    
    % Create an interpolated model function
    % This returns a new function that is a weighted combination of the two closest models
    model_fun = @(x, y) (1-w) * model1(x, y) + w * model2(x, y);
end