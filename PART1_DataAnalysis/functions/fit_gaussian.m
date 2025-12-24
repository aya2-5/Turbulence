function G=fit_gaussian(du, pdf_du)
    % du     : vector of velocity increments [m]
    % pdf_du : probability density function (PDF)
    % G      : A gaussian distribution fit G(du)
     % Estimate mean (μ) and variance (σ^2) from the PDF
    mu    = trapz(du, du .* pdf_du);
    sigma2 = trapz(du, (du - mu).^2 .* pdf_du);
    sigma = sqrt(sigma2);
    
    % Build Gaussian with those parameters
    G = (1 ./ (sigma * sqrt(2*pi))) .* exp( - (du - mu).^2 ./ (2 * sigma2) );
end
    