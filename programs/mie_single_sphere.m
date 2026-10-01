function result = mie_single_sphere(theta,diameter,lambda0,nParticle,nMedium)
%MIE_SINGLE_SPHERE Far-field Mie scattering; all lengths in meters, angles rad.
% lambda0 is VACUUM wavelength; nParticle is ABSOLUTE n-i*kappa.
% S1 perpendicular / S2 parallel to the scattering plane, raw series scale.
% differentialCrossSection has units m^2/sr for unpolarized incidence.
% PhaseFunction integrates to one over 4*pi when scattering is nonzero.
    if nargin<5, nMedium=1; end
    validateattributes(theta,{'numeric'},{'vector','nonempty','real','finite','>=',0,'<=',pi});
    validateattributes(diameter,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(lambda0,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(nMedium,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(nParticle,{'numeric'},{'scalar','finite'});
    theta = theta(:);
    k = 2*pi*nMedium/lambda0;
    x = k*diameter/2;
    m = nParticle/nMedium;
    [a,b,info] = mie_coefficients(m,x);
    mu = cos(theta);
    previousPi = zeros(size(mu));
    currentPi = ones(size(mu));
    S1 = complex(zeros(size(mu)));
    S2 = S1;
    for n = 1:numel(a)
        tau = n*mu.*currentPi-(n+1)*previousPi;
        factor = (2*n+1)/(n*(n+1));
        S1 = S1+factor*(a(n)*currentPi+b(n)*tau);
        S2 = S2+factor*(a(n)*tau+b(n)*currentPi);
        nextPi = ((2*n+1)*mu.*currentPi-(n+1)*previousPi)/n;
        previousPi = currentPi;
        currentPi = nextPi;
    end
    n = (1:numel(a)).';
    if x==0
        Qext=0; Qsca=0; Qback=0; g=0;
    else
        Qext = 2/x^2*sum((2*n+1).*real(a+b));
        Qsca = 2/x^2*sum((2*n+1).*(abs(a).^2+abs(b).^2));
        Qback = abs(sum((2*n+1).*(-1).^n.*(a-b)))^2/x^2;
        g=0;
        if Qsca>0
            adjacent = n(1:end-1).*(n(1:end-1)+2)./(n(1:end-1)+1);
            sumAdjacent = sum(adjacent.*real(a(1:end-1).*conj(a(2:end)) ...
                +b(1:end-1).*conj(b(2:end))));
            sumCross = sum((2*n+1)./(n.*(n+1)).*real(a.*conj(b)));
            g = 4*(sumAdjacent+sumCross)/(x^2*Qsca);
        end
    end
    Qabs = Qext-Qsca;
    if Qabs < -1e-10*max(1,Qext)
        error('Mie:NegativeAbsorption','Passive sphere produced negative absorption.');
    end
    area = pi*diameter^2/4;
    iPerpendicular = abs(S1).^2;
    iParallel = abs(S2).^2;
    intensityFunction = (iPerpendicular+iParallel)/2;
    dcs = intensityFunction/k^2;
    Csca = area*Qsca;
    phase = zeros(size(dcs));
    if Csca>0, phase=dcs/Csca; end
    result = struct('theta',theta,'diameter',diameter,'lambda0',lambda0, ...
        'nParticle',nParticle,'nMedium',nMedium,'m',m,'x',x,'k',k, ...
        'a',a,'b',b,'S1',S1,'S2',S2,'Qext',Qext,'Qsca',Qsca, ...
        'Qabs',Qabs,'Qback',Qback,'g',g,'Cext',area*Qext, ...
        'Csca',Csca,'Cabs',area*Qabs,'iPerpendicular',iPerpendicular, ...
        'iParallel',iParallel,'intensityFunction',intensityFunction, ...
        'differentialCrossSection',dcs,'phaseFunction',phase,'info',info);
end
