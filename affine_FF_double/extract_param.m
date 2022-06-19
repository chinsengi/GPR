function [sigmaE, lE, snE,...
              sigmaI, lI, snI, innoise]...
                                     = extract_param(X)
    % it's a hack, so that we don't need to know if there is input noise
    X = [X; 0];
    X = exp(X);
    [sigmaE, lE, snE, sigmaI, lI, snI, innoise] =...
                                              feval(@(x) x{:}, num2cell(X));
end