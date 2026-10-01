function [a,b,info] = mie_coefficients(m,x)
%MIE_COEFFICIENTS External coefficients of a homogeneous nonmagnetic sphere.
% Convention: exp(+i*w*t), outgoing h_n^(2), passive m = n - i*kappa.
% m is relative to a real, nonabsorbing host; x = pi*nMedium*d/lambda0.
% Supported x: zero, or 1e-6 <= x <= 1e4. No perfect-conductor sentinel.
% Formulas are evaluated using real-argument Bessel functions and a continued
% fraction for D_n(mx); complex Bessel functions at mx are not evaluated.

    validateattributes(m,{'numeric'},{'scalar','finite'},mfilename,'m');
    validateattributes(x,{'numeric'},{'scalar','real','finite','nonnegative'},mfilename,'x');
    if real(m)<=0 || imag(m)>0
        error('Mie:IndexConvention','Require real(m)>0 and imag(m)<=0 (n-i*kappa).');
    end
    if x>1e4 || (x>0 && x<1e-6)
        error('Mie:SizeRange','Require x=0 or 1e-6 <= x <= 1e4.');
    end
    nmax = max(2,ceil(x+4.05*x^(1/3)+2));
    a = complex(zeros(nmax,1));
    b = a;
    info = struct('nTerms',nmax,'continuedFractionIterations',0, ...
        'timeConvention','exp(+i*omega*t)','outgoingHankelKind',2);
    if x==0 || m==1
        return;
    end

    z = m*x;
    [dn,iterations] = log_derivative_cf(z,nmax);
    info.continuedFractionIterations = iterations;
    D = complex(zeros(nmax,1));
    D(nmax) = dn;
    for n = nmax:-1:2
        D(n-1) = n/z - 1/(D(n)+n/z);
    end
    orders = (0:nmax).';
    psi = sqrt(pi*x/2)*besselj(orders+0.5,x);
    chi = -sqrt(pi*x/2)*bessely(orders+0.5,x);
    xi = psi+1i*chi;  % x*h_n^(2)(x) = x*j_n(x) - i*x*y_n(x)
    n = (1:nmax).';
    ta = D/m+n/x;
    tb = m*D+n/x;
    a = (ta.*psi(2:end)-psi(1:end-1))./(ta.*xi(2:end)-xi(1:end-1));
    b = (tb.*psi(2:end)-psi(1:end-1))./(tb.*xi(2:end)-xi(1:end-1));
    if any(~isfinite(a)) || any(~isfinite(b))
        error('Mie:NonfiniteCoefficients','Coefficient calculation failed for this m,x.');
    end
end

function [dn,iterations] = log_derivative_cf(z,n)
% D_n = J_(n-1/2)/J_(n+1/2) - n/z.
% The Bessel ratio is (2n+1)/z - 1/((2n+3)/z - 1/(...)).
% Modified Lentz iteration; fail explicitly if the fraction does not converge.
    tiny = 1e-300;
    f = (2*n+1)/z;
    if abs(f)<tiny, f=tiny; end
    C = f;
    inverseDenominator = 0;
    for iterations = 1:100000
        bj = (2*(n+iterations)+1)/z;
        denominator = bj-inverseDenominator;
        if abs(denominator)<tiny, denominator=tiny; end
        C = bj-1/C;
        if abs(C)<tiny, C=tiny; end
        inverseDenominator = 1/denominator;
        delta = C*inverseDenominator;
        f = f*delta;
        if abs(delta-1)<1e-13
            dn = f-n/z;
            return;
        end
    end
    error('Mie:NoConvergence','Logarithmic derivative fraction did not converge.');
end
