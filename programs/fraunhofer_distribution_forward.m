function result = fraunhofer_distribution_forward(theta,diameter,psd,lambda0,nMedium,options,variant,cached)
%FRAUNHOFER_DISTRIBUTION_FORWARD Same PSD weighting/units as the Mie forward model.
% Optional options fields: DistributionBasis, DistributionSampling,
% IncidentIrradiance, Distance, TotalNumber. Scalar approximation has no azimuth.
% Polarization/Azimuth from Mie options are accepted as metadata only.
% Optional cached: reuse a previous result with identical optical/grid settings.
    if nargin<5 || isempty(nMedium), nMedium=1; end
    if nargin<6 || isempty(options), options=struct(); end
    if nargin<7 || isempty(variant), variant='sphere'; end
    defaults=struct('DistributionBasis','number','DistributionSampling','density', ...
        'IncidentIrradiance',1,'Distance',1,'TotalNumber',[], ...
        'Polarization','unpolarized','Azimuth',0);
    if ~isstruct(options) || ~isscalar(options), error('Mie:Options','options must be a scalar struct.'); end
    fields=fieldnames(options);
    for j=1:numel(fields)
        if ~isfield(defaults,fields{j}), error('Mie:UnknownOption','Unknown option: %s.',fields{j}); end
        defaults.(fields{j})=options.(fields{j});
    end
    options=defaults;
    variant=validatestring(variant,{'sphere','paraxial'});
    validateattributes(theta,{'numeric'},{'vector','nonempty','real','finite','>=',0,'<',pi/2});
    validateattributes(lambda0,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(nMedium,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(options.IncidentIrradiance,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(options.Distance,{'numeric'},{'scalar','real','finite','positive'});
    if ~isempty(options.TotalNumber)
        validateattributes(options.TotalNumber,{'numeric'},{'scalar','real','finite','nonnegative'});
    end
    [fractions,weights,scale]=particle_distribution_weights(diameter,psd, ...
        options.DistributionBasis,options.DistributionSampling);
    cacheKey=struct('version','fraunhofer-kernel-v1','theta',theta(:),'diameter',diameter(:), ...
        'lambda0',lambda0,'nMedium',nMedium,'variant',variant);
    cacheHit=nargin>=8 && ~isempty(cached);
    if cacheHit
        if ~isstruct(cached) || ~isscalar(cached) || ~all(isfield(cached,{'cacheKey','kernel','info'})) || ...
                ~isequaln(cached.cacheKey,cacheKey)
            error('Fraunhofer:CacheMismatch','Cached optical settings do not match this calculation.');
        end
        kernel=cached.kernel; info=cached.info;
        if ~isequal(size(kernel),[numel(theta),numel(diameter)])
            error('Fraunhofer:CacheMismatch','Cached matrix size does not match this calculation.');
        end
    else
        [kernel,info]=fraunhofer_scattering_kernel(theta,diameter,lambda0,nMedium,variant);
    end
    dcs=kernel*fractions; peak=max(dcs); relative=zeros(size(dcs));
    if peak>0, relative=dcs/peak; end
    irradiance=[];
    if ~isempty(options.TotalNumber)
        irradiance=options.IncidentIrradiance*options.TotalNumber*dcs/options.Distance^2;
    end
    result=struct('theta',theta(:),'diameter',diameter(:),'psd',psd(:), ...
        'options',options,'lambda0',lambda0,'nMedium',nMedium,'numberFractions',fractions, ...
        'numberScale',scale,'kernel',kernel,'distributionKernel',kernel.*weights.', ...
        'differentialCrossSection',dcs,'relativeIntensity',relative,'irradiance',irradiance, ...
        'meanNominalDiffractionCrossSection',info.nominalDiffractionCrossSection*fractions, ...
        'info',info,'cacheKey',cacheKey,'cacheHit',cacheHit);
end
