function report = test_mie_scattering()
%TEST_MIE_SCATTERING Independent reference and physical consistency checks.
    folder=fileparts(fileparts(mfilename('fullpath')));
    data=jsondecode(fileread(fullfile(folder,'tests','fixtures','miepython_reference.json')));
    theta=data.theta(:);
    errors=zeros(numel(data.cases),3);
    for j=1:numel(data.cases)
        c=data.cases(j);
        m=complex(c.mReal,c.mImag);
        % lambda0=pi makes diameter=x for a host with n=1.
        s=mie_single_sphere(theta,c.x,pi,m,1);
        actual=[s.Qext;s.Qsca;s.Qback;s.g];
        expected=c.efficiencies(:);
        % g tends to zero as x^2: compare near-zero g in absolute accuracy,
        % while retaining a relative comparison for nonzero efficiencies.
        denominator=max(abs(expected),[1e-14;1e-14;1e-14;1e-8]);
        errors(j,1)=max(abs(actual-expected)./denominator);
        expectedS1=complex(c.S1Real(:),c.S1Imag(:));
        expectedS2=complex(c.S2Real(:),c.S2Imag(:));
        scale=max(max(abs([expectedS1;expectedS2])),1e-30);
        if m==1, scale=1; end  % exact zero; reference leaves roundoff residuals
        errors(j,2)=max(abs([s.S1-expectedS1;s.S2-expectedS2]))/scale;
        % Public reference coefficients are conjugated relative to its S1_S2.
        expectedA=conj(complex(c.aReal(:),c.aImag(:)));
        expectedB=conj(complex(c.bReal(:),c.bImag(:)));
        scale=max(max(abs([expectedA;expectedB])),1e-30);
        if m==1, scale=1; end
        errors(j,3)=max(abs([s.a(1:2)-expectedA;s.b(1:2)-expectedB]))/scale;
        assert(abs(s.S1(1)-s.S2(1))<1e-10*max(1,abs(s.S1(1))));
        assert(abs(s.S1(end)+s.S2(end))<1e-10*max(1,abs(s.S1(end))));
        assert(s.Qsca>=0 && s.Qabs>=-1e-10*max(1,s.Qext));
    end
    % Reference uses a small-sphere approximation for efficiencies below |mx|=.1;
    % MATLAB evaluates the full series. Compare those with an explicit tolerance.
    assert(max(errors(:,1))<2e-5,'Efficiency reference comparison failed.');
    assert(max(errors(:,2))<1e-6,'Complex amplitude reference comparison failed.');
    assert(max(errors(:,3))<1e-6,'Complex coefficient reference comparison failed.');

    % Gauss-Legendre integration: polynomials from the truncated angular series
    % integrate exactly at sufficient order; no forward-peak sampling ambiguity.
    physicalCases=[1.5,0,1;1.5,-0.1,10;0.75,0,10;1.5,-0.1,100; ...
        1.591,0,pi*200e-6/450e-9];
    energyError=zeros(size(physicalCases,1),1);
    for j=1:size(physicalCases,1)
        m=complex(physicalCases(j,1),physicalCases(j,2)); x=physicalCases(j,3);
        [a,~,~]=mie_coefficients(m,x);
        [mu,w]=gauss_legendre(numel(a)+2);
        s=mie_single_sphere(acos(mu),x,pi,m,1);
        angular=2*pi*sum(w.*s.intensityFunction);
        expected=pi*x^2*s.Qsca;
        energyError(j)=abs(angular-expected)/expected;
        assert(energyError(j)<1e-10,'Scattering integral failed.');
        forward=mie_single_sphere(0,x,pi,m,1);
        assert(abs(4*real(forward.S1)/x^2-s.Qext)<1e-10);
        assert(abs(2*pi*sum(w.*s.phaseFunction)-1)<1e-10);
        if imag(m)==0, assert(abs(s.Qabs)<1e-10); end
    end
    ray=mie_single_sphere([0;pi/2;pi],1e-4,pi,1.5,1);
    rayExpected=8/3*ray.x^4*abs((1.5^2-1)/(1.5^2+2))^2;
    assert(abs(ray.Qsca/rayExpected-1)<1e-7);
    assert(abs(ray.iParallel(2)/ray.iPerpendicular(2))<1e-12);
    assert(max(abs(ray.intensityFunction/ray.intensityFunction(2)-[2;1;2]))<1e-7);
    zero=mie_single_sphere([0;pi],0,633e-9,1.5,1);
    assert(all(zero.intensityFunction==0) && zero.Qsca==0);

    e=data.ensemble;
    options=struct('DistributionBasis','volume','DistributionSampling','density');
    mixture=mie_distribution_forward(theta,e.diameter,e.volumeDensity,e.lambda0, ...
        complex(e.nParticleReal,e.nParticleImag),e.nMedium,options);
    ensembleError=max(abs(mixture.differentialCrossSection-e.differentialCrossSection(:))) ...
        /max(e.differentialCrossSection);
    assert(ensembleError<1e-6,'Ensemble reference comparison failed.');
    d=e.diameter(:); volume=pi*d.^3/6;
    numberOptions=struct('DistributionBasis','number');
    equivalent=mie_distribution_forward(theta,d,e.volumeDensity(:)./volume,e.lambda0, ...
        complex(e.nParticleReal,e.nParticleImag),e.nMedium,numberOptions);
    assert(norm(equivalent.differentialCrossSection-mixture.differentialCrossSection) ...
        <1e-12*norm(mixture.differentialCrossSection));
    assert(norm(mixture.distributionKernel*mixture.psd/mixture.numberScale ...
        -mixture.differentialCrossSection)<1e-12*norm(mixture.differentialCrossSection));
    binOptions=struct('DistributionSampling','bin');
    bin=mie_distribution_forward(theta,d,mixture.numberFractions,e.lambda0, ...
        complex(e.nParticleReal,e.nParticleImag),e.nMedium,binOptions);
    assert(norm(bin.differentialCrossSection-mixture.differentialCrossSection) ...
        <1e-12*norm(mixture.differentialCrossSection));

    single=mie_single_sphere(theta,2e-6,633e-9,1.59-0.01i,1.33);
    mono=mie_distribution_forward(theta,2e-6,1,633e-9,1.59-0.01i,1.33,binOptions);
    assert(norm(single.differentialCrossSection-mono.differentialCrossSection)==0);
    options=struct('DistributionSampling','bin','Polarization','linear','Azimuth',0, ...
        'TotalNumber',200,'Distance',2,'IncidentIrradiance',3);
    linear=mie_distribution_forward(theta,2e-6,1,633e-9,1.59-0.01i,1.33,options);
    assert(norm(linear.differentialCrossSection-single.iParallel/single.k^2)==0);
    assert(norm(linear.irradiance-150*linear.differentialCrossSection)==0);
    options.Azimuth=pi/2;
    linear=mie_distribution_forward(theta,2e-6,1,633e-9,1.59-0.01i,1.33,options);
    assert(norm(linear.differentialCrossSection-single.iPerpendicular/single.k^2) ...
        <1e-12*norm(single.differentialCrossSection));

    % Invalid physics/data must fail visibly rather than be silently coerced.
    must_error(@()mie_coefficients(1.5+0.1i,1),'Mie:IndexConvention');
    must_error(@()mie_coefficients(1.5,1e-7),'Mie:SizeRange');
    must_error(@()mie_distribution_forward(theta,[1;2]*1e-6,[0;0],633e-9,1.5,1), ...
        'Mie:ZeroPSD');
    must_error(@()mie_distribution_forward(theta,1e-6,1,633e-9,1.5,1), ...
        'Mie:DensityGrid');
    must_error(@()mie_distribution_forward(theta,[2;1]*1e-6,[1;1],633e-9,1.5,1), ...
        'Mie:Grid');
    must_error(@()mie_distribution_forward(theta,[1;2]*1e-6,[1;1],633e-9,1.5,1, ...
        struct('Unknown',1)),'Mie:UnknownOption');

    report=struct('caseCount',numel(data.cases),'referenceVersion',data.packageVersion, ...
        'referenceCommit',data.commit,'maxEfficiencyScaledError',max(errors(:,1)), ...
        'maxAmplitudeScaledError',max(errors(:,2)), ...
        'maxCoefficientScaledError',max(errors(:,3)), ...
        'maxEnergyRelativeError',max(energyError),'ensembleScaledError',ensembleError);
    disp(report);
    fprintf('MIE_VALIDATION_OK\n');
end

function [mu,w]=gauss_legendre(n)
    j=(1:n-1).'; beta=j./sqrt(4*j.^2-1);
    [V,L]=eig(diag(beta,1)+diag(beta,-1));
    [mu,index]=sort(diag(L)); w=2*(V(1,index).').^2;
end

function must_error(action,identifier)
    try
        action();
    catch exception
        assert(strcmp(exception.identifier,identifier),'Unexpected error identifier.');
        return;
    end
    error('Mie:ExpectedError','Expected error %s was not raised.',identifier);
end
