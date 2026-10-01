function report = test_scattering_comparison()
%TEST_SCATTERING_COMPARISON Verify distributions, diffraction scale and cache.
    d=linspace(10,200,5001).'*1e-6;
    definitions={struct('Scale',80e-6,'Shape',3), ...
        struct('Mean',70e-6,'Std',20e-6),struct('Median',70e-6,'LogSigma',0.3), ...
        struct('Medians',[40,120]*1e-6,'LogSigmas',[0.16,0.18],'VolumeFractions',[0.3,0.7])};
    types={'rr','normal','lognormal','bimodal'}; cdfErrors=zeros(4,1);
    for j=1:4
        g=generate_particle_distribution(d,types{j},definitions{j});
        assert(all(isfinite(g.density)) && all(g.density>=0));
        assert(abs(trapz(d,g.density)-1)<1e-12 && all(diff(g.cdf)>=0));
        switch j
            case 1, cdf=1-exp(-(d/80e-6).^3); cdf=truncate_cdf(cdf);
            case 2, cdf=0.5*(1+erf((d-70e-6)/(20e-6*sqrt(2)))); cdf=truncate_cdf(cdf);
            case 3, cdf=0.5*(1+erf(log(d/70e-6)/(0.3*sqrt(2)))); cdf=truncate_cdf(cdf);
            case 4
                c1=truncate_cdf(0.5*(1+erf(log(d/40e-6)/(0.16*sqrt(2)))));
                c2=truncate_cdf(0.5*(1+erf(log(d/120e-6)/(0.18*sqrt(2)))));
                cdf=0.3*c1+0.7*c2;
        end
        cdfErrors(j)=max(abs(g.cdf-cdf));
        assert(cdfErrors(j)<2e-6,'Distribution differs from its analytic truncated CDF.');
    end
    % A mode is not the R-R scale or the lognormal median. Check the new input
    % route against analytic CDFs and the actual locations of density maxima.
    rrMode=generate_particle_distribution(d,'rr',struct('Mode',100e-6,'Shape',3));
    logMode=generate_particle_distribution(d,'lognormal',struct('Mode',100e-6,'LogSigma',0.3));
    [~,rrIndex]=max(rrMode.density); [~,logIndex]=max(logMode.density);
    modeErrors=abs([d(rrIndex);d(logIndex)]-100e-6);
    assert(max(modeErrors)<=max(diff(d)));
    scale=100e-6/(2/3)^(1/3);
    rrCdf=truncate_cdf(1-exp(-(d/scale).^3));
    logCdf=truncate_cdf(0.5*(1+erf(log(d/(100e-6*exp(0.09)))/(0.3*sqrt(2)))));
    assert(max(abs(rrMode.cdf-rrCdf))<2e-6 && max(abs(logMode.cdf-logCdf))<2e-6);
    twoModes=generate_particle_distribution(d,'bimodal', ...
        struct('Modes',[80,130]*1e-6,'LogSigmas',[0.15,0.10],'VolumeFractions',[0.5,0.5]));
    [~,indices]=max(twoModes.componentDensities,[],1);
    assert(max(abs(d(indices)-[80;130]*1e-6))<=max(diff(d)));
    must_error(@()generate_particle_distribution(d,'rr',struct('Mode',100e-6,'Shape',1)),'PSD:NoInteriorMode');
    lambda0=450e-9; diameter=100e-6; k=2*pi/lambda0; a=diameter/2; x=k*a;
    theta=[0;0.02;0.1;0.2;1]*pi/180;
    [K,~]=fraunhofer_scattering_kernel(theta,diameter,lambda0,1,'sphere');
    assert(abs(K(1)/(k^2*a^4/4)-1)<1e-14);
    % Independent integration over a circular aperture using J0.
    integralDcs=zeros(size(theta));
    for j=1:numel(theta)
        radial=integral(@(r)r.*besselj(0,k*sin(theta(j))*r),0,a,'AbsTol',1e-24,'RelTol',1e-12);
        amplitude=k*radial*(1+cos(theta(j)))/2;
        integralDcs(j)=abs(amplitude)^2;
    end
    apertureError=max(abs(K-integralDcs))/max(K);
    assert(apertureError<1e-10);
    firstZero=asin(3.8317059702075125/x);
    zeroKernel=fraunhofer_scattering_kernel(firstZero,diameter,lambda0,1);
    assert(zeroKernel/K(1)<1e-25);
    % First Airy side lobe is about 1.75% of the main peak; it must be present
    % without artificially adding peaks or smoothing the computed intensity.
    sideAngle=fminbnd(@(t)-fraunhofer_scattering_kernel(t,diameter,lambda0,1), ...
        firstZero,asin(7.015586669815619/x),optimset('TolX',1e-12));
    firstSideFraction=fraunhofer_scattering_kernel(sideAngle,diameter,lambda0,1)/K(1);
    assert(firstSideFraction>0.0174 && firstSideFraction<0.0176);
    doubled=fraunhofer_scattering_kernel(0,2*diameter,lambda0,1);
    assert(abs(doubled/K(1)-16)<1e-12); % forward scale is d^4
    red=fraunhofer_scattering_kernel(0,diameter,632e-9,1);
    assert(abs(red/K(1)-(450/632)^2)<1e-14);
    paraxial=fraunhofer_scattering_kernel(theta,diameter,lambda0,1,'paraxial');
    sphere=fraunhofer_scattering_kernel(theta,diameter,lambda0,1,'sphere');
    assert(abs(paraxial(1)-sphere(1))==0);
    % Matched index produces zero Mie scattering: size alone cannot guarantee
    % that a refractive sphere is represented by an opaque-disk approximation.
    matched=mie_single_sphere(theta,diameter,lambda0,1,1);
    assert(all(matched.differentialCrossSection==0) && K(1)>0);

    smallD=linspace(10,200,31).'*1e-6; angles=linspace(0,2,121).'*pi/180;
    p1=generate_particle_distribution(smallD,'normal',definitions{2});
    p2=generate_particle_distribution(smallD,'bimodal',definitions{4});
    opt=struct('DistributionBasis','volume');
    base=mie_distribution_forward(angles,smallD,p1.density,lambda0,1.591,1,opt);
    reused=mie_distribution_forward(angles,smallD,p2.density,lambda0,1.591,1,opt,base);
    direct=mie_distribution_forward(angles,smallD,p2.density,lambda0,1.591,1,opt);
    cacheError=norm(reused.differentialCrossSection-direct.differentialCrossSection) ...
        /norm(direct.differentialCrossSection);
    assert(reused.cacheHit && ~direct.cacheHit && cacheError<1e-14);
    f=fraunhofer_distribution_forward(angles,smallD,p2.density,lambda0,1,opt);
    assert(isequal(f.numberFractions,reused.numberFractions));
    assert(norm(f.distributionKernel*p2.density/f.numberScale-f.differentialCrossSection) ...
        <1e-12*norm(f.differentialCrossSection));
    must_error(@()mie_distribution_forward(angles,smallD,p1.density,632e-9,1.591,1,opt,base),'Mie:CacheMismatch');
    must_error(@()fraunhofer_scattering_kernel(pi,diameter,lambda0,1),'');

    % Known factor-of-two curves distinguish amplitude differences from shape.
    signal=exp(-angles); metrics=compare_scattering_models(angles,signal,2*signal,pi/180);
    assert(abs(metrics.angularL2Percent-100)<1e-10);
    assert(metrics.independentPeakShapeL2Percent<1e-10);
    assert(abs(metrics.fullConeDifferencePercent-100)<1e-10);
    % Check the main workflow with four PSDs, two wavelengths and no file I/O.
    overrides=struct('diameterPoints',31,'anglePoints',121,'makePlots',false,'saveData',false,'saveFigures',false, ...
        'showPeakDiagnostics',false);
    s=main_mie_scattering(overrides);
    assert(height(s.metrics)==8 && numel(s.distributions)==4);
    assert(~s.mieResults{1,1}.cacheHit && all(cellfun(@(r)r.cacheHit,s.mieResults(2:end,1))));
    report=struct('maxTruncatedCdfError',max(cdfErrors),'apertureScaledError',apertureError, ...
        'cacheRelativeError',cacheError,'mainComparisonCount',height(s.metrics), ...
        'maxModeLocationError_um',max(modeErrors)*1e6,'firstDiffractionSidePeakFraction',firstSideFraction);
    disp(report); fprintf('SCATTERING_COMPARISON_VALIDATION_OK\n');
end

function cdf=truncate_cdf(cdf)
    cdf=(cdf-cdf(1))/(cdf(end)-cdf(1));
end

function must_error(action,identifier)
    try
        action();
    catch exception
        if ~isempty(identifier), assert(strcmp(exception.identifier,identifier)); end
        return;
    end
    error('Test:ExpectedError','Expected an error was not raised.');
end
