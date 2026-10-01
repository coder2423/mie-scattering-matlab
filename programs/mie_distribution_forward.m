function result = mie_distribution_forward(theta,diameter,psd,lambda0,nParticle,nMedium,options,cached)
%MIE_DISTRIBUTION_FORWARD Incoherent single-scattering average of a PSD.
% theta rad; diameter/lambda0 m; nParticle absolute n-i*kappa; host real.
% DistributionBasis: 'number' (default) / 'volume'.
% DistributionSampling: 'density' per meter of diameter (default) / 'bin'
% (probability/fraction already integrated over bins; no width multiplier).
% DENSITY nodes are strictly increasing and use trapezoidal quadrature;
% one-point monodisperse input must use BIN sampling.
% Polarization: 'unpolarized' (default), 'parallel', 'perpendicular', 'linear'.
% Azimuth rad: for linear incidence, angle from incident E to scattering plane.
% TotalNumber: optional total illuminated particle count (not concentration).
% Irradiance at Distance is only returned when TotalNumber is provided.
% Each particle uses the same incident irradiance; negligible attenuation.
% No per-diameter phase-function normalization is used in the mixture.
% Optional cached: previous result with identical optical/grid/polarization
% settings; reuse its single-particle kernel when only the PSD changes.
    if nargin<6 || isempty(nMedium), nMedium=1; end
    if nargin<7 || isempty(options), options=struct(); end
    defaults = struct('DistributionBasis','number', ...
        'DistributionSampling','density','Polarization','unpolarized', ...
        'Azimuth',0,'IncidentIrradiance',1,'Distance',1,'TotalNumber',[]);
    if ~isstruct(options) || ~isscalar(options)
        error('Mie:Options','options must be a scalar struct.');
    end
    names = fieldnames(options);
    for j=1:numel(names)
        if ~isfield(defaults,names{j})
            error('Mie:UnknownOption','Unknown option: %s.',names{j});
        end
        defaults.(names{j})=options.(names{j});
    end
    options=defaults;
    basis=validatestring(options.DistributionBasis,{'number','volume'});
    sampling=validatestring(options.DistributionSampling,{'density','bin'});
    polarization=validatestring(options.Polarization, ...
        {'unpolarized','parallel','perpendicular','linear'});
    validateattributes(theta,{'numeric'},{'vector','nonempty','real','finite','>=',0,'<=',pi});
    validateattributes(diameter,{'numeric'},{'vector','nonempty','real','finite','positive'});
    validateattributes(psd,{'numeric'},{'vector','nonempty','real','finite','nonnegative'});
    validateattributes(options.IncidentIrradiance,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(options.Distance,{'numeric'},{'scalar','real','finite','positive'});
    if ~isempty(options.TotalNumber)
        validateattributes(options.TotalNumber,{'numeric'},{'scalar','real','finite','nonnegative'});
    end
    diameter=diameter(:); psd=psd(:); theta=theta(:);
    [numberFractions,weights,numberScale]=particle_distribution_weights(diameter,psd,basis,sampling);
    validateattributes(lambda0,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(nMedium,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(nParticle,{'numeric'},{'scalar','finite'});
    if real(nParticle)<=0 || imag(nParticle)>0
        error('Mie:IndexConvention','Require real(nParticle)>0 and imag(nParticle)<=0.');
    end
    phi=options.Azimuth;
    validateattributes(phi,{'numeric'},{'real','finite','nonempty'});
    if ~isscalar(phi) && (~isvector(phi) || numel(phi)~=numel(theta))
        error('Mie:Azimuth','Azimuth must be scalar or have the same length as theta.');
    end
    phi=phi(:);
    cacheKey=struct('version','mie-kernel-v1','theta',theta,'diameter',diameter, ...
        'lambda0',lambda0,'nParticle',nParticle,'nMedium',nMedium, ...
        'polarization',polarization,'azimuth',phi);
    cacheHit=nargin>=8 && ~isempty(cached);
    if cacheHit
        if ~isstruct(cached) || ~isscalar(cached) || ...
                ~all(isfield(cached,{'cacheKey','kernel','efficiencies','crossSections'})) || ~isequaln(cached.cacheKey,cacheKey)
            error('Mie:CacheMismatch','Cached optical settings do not match this calculation.');
        end
        kernel=cached.kernel; efficiency=cached.efficiencies; crossSections=cached.crossSections;
        if ~isequal(size(kernel),[numel(theta),numel(diameter)]) || ...
                ~isequal(size(efficiency),[numel(diameter),4]) || ...
                ~isequal(size(crossSections),[numel(diameter),3])
            error('Mie:CacheMismatch','Cached matrix sizes do not match this calculation.');
        end
    else
    kernel=zeros(numel(theta),numel(diameter));
    efficiency=zeros(numel(diameter),4);
    crossSections=zeros(numel(diameter),3);
    for j=1:numel(diameter)
        sphere=mie_single_sphere(theta,diameter(j),lambda0,nParticle,nMedium);
        switch polarization
            case 'unpolarized', intensity=sphere.intensityFunction;
            case 'parallel', intensity=sphere.iParallel;
            case 'perpendicular', intensity=sphere.iPerpendicular;
            case 'linear'
                intensity=sphere.iParallel.*cos(phi).^2 ...
                    +sphere.iPerpendicular.*sin(phi).^2;
        end
        kernel(:,j)=intensity/sphere.k^2;
        efficiency(j,:)=[sphere.Qext,sphere.Qsca,sphere.Qabs,sphere.g];
        crossSections(j,:)=[sphere.Cext,sphere.Csca,sphere.Cabs];
    end
    end
    dcs=kernel*numberFractions;
    meanCrossSections=crossSections.'*numberFractions;
    peak=max(dcs);
    relative=zeros(size(dcs));
    if peak>0, relative=dcs/peak; end
    % For density samples K includes quadrature. Volume samples also include 1/V.
    % K*psd / numberScale equals the mean per-particle differential cross section.
    distributionKernel=kernel.*weights.';
    irradiance=[];
    if ~isempty(options.TotalNumber)
        irradiance=options.IncidentIrradiance*options.TotalNumber*dcs/options.Distance^2;
    end
    phase=[];
    if strcmp(polarization,'unpolarized')
        phase=zeros(size(dcs));
        if meanCrossSections(2)>0, phase=dcs/meanCrossSections(2); end
    end
    result=struct('theta',theta,'diameter',diameter,'psd',psd, ...
        'options',options,'lambda0',lambda0,'nParticle',nParticle, ...
        'nMedium',nMedium,'numberFractions',numberFractions, ...
        'numberScale',numberScale,'kernel',kernel, ...
        'distributionKernel',distributionKernel,'efficiencies',efficiency, ...
        'crossSections',crossSections,'cacheKey',cacheKey,'cacheHit',cacheHit, ...
        'meanCext',meanCrossSections(1),'meanCsca',meanCrossSections(2), ...
        'meanCabs',meanCrossSections(3),'differentialCrossSection',dcs, ...
        'relativeIntensity',relative,'phaseFunction',phase,'irradiance',irradiance);
end
