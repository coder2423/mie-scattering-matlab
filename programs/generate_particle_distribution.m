function distribution = generate_particle_distribution(diameter,type,parameters)
%GENERATE_PARTICLE_DISTRIBUTION Truncated volume/number density on a D grid.
% All size parameters are meters; LogSigma/Shape/VolumeFractions dimensionless.
% The caller declares the basis. Main uses all four as VOLUME densities.
% rr: Mode (or Scale),Shape; normal: Mean,Std;
% lognormal: Mode (or Median),LogSigma;
% bimodal: Modes (or Medians),LogSigmas,VolumeFractions.
% Mixture fractions refer to fractions WITHIN the configured diameter range.
    validateattributes(diameter,{'numeric'},{'vector','real','finite','positive'});
    diameter=diameter(:);
    if numel(diameter)<2 || any(diff(diameter)<=0)
        error('PSD:Grid','Density grid needs at least two strictly increasing nodes.');
    end
    type=validatestring(type,{'rr','normal','lognormal','bimodal'});
    if ~isstruct(parameters) || ~isscalar(parameters)
        error('PSD:Parameters','parameters must be a scalar struct.');
    end
    componentDensities=[]; fractions=[]; effectiveParameters=parameters;
    switch type
        case 'rr'
            shape=positive_parameter(parameters,'Shape');
            if isfield(parameters,'Mode')
                if isfield(parameters,'Scale'), error('PSD:AmbiguousParameter','Set Mode or Scale, not both.'); end
                mode=positive_parameter(parameters,'Mode');
                if shape<=1, error('PSD:NoInteriorMode','R-R Mode requires Shape>1.'); end
                scale=mode/((shape-1)/shape)^(1/shape);
                effectiveParameters.Scale=scale;
            else
                scale=positive_parameter(parameters,'Scale');
            end
            logRatio=log(diameter/scale);
            logDensity=log(shape/scale)+(shape-1)*logRatio-exp(shape*logRatio);
            density=normalize_log_density(diameter,logDensity);
            name='Rosin-Rammler';
        case 'normal'
            center=positive_parameter(parameters,'Mean');
            sigma=positive_parameter(parameters,'Std');
            density=normalize_log_density(diameter,-0.5*((diameter-center)/sigma).^2);
            name='Normal';
        case 'lognormal'
            sigma=positive_parameter(parameters,'LogSigma');
            if isfield(parameters,'Mode')
                if isfield(parameters,'Median'), error('PSD:AmbiguousParameter','Set Mode or Median, not both.'); end
                center=positive_parameter(parameters,'Mode')*exp(sigma^2);
                effectiveParameters.Median=center;
            else
                center=positive_parameter(parameters,'Median');
            end
            density=normalize_log_density(diameter,-log(diameter) ...
                -0.5*(log(diameter/center)/sigma).^2);
            name='Lognormal';
        case 'bimodal'
            sigmas=parameters.LogSigmas(:);
            validateattributes(sigmas,{'numeric'},{'real','finite','positive','numel',2});
            if isfield(parameters,'Modes')
                if isfield(parameters,'Medians'), error('PSD:AmbiguousParameter','Set Modes or Medians, not both.'); end
                modes=parameters.Modes(:);
                validateattributes(modes,{'numeric'},{'real','finite','positive','numel',2});
                centers=modes.*exp(sigmas.^2);
                effectiveParameters.Medians=centers;
            else
                centers=parameters.Medians(:);
            end
            fractions=parameters.VolumeFractions(:);
            validateattributes(centers,{'numeric'},{'real','finite','positive','numel',2});
            validateattributes(sigmas,{'numeric'},{'real','finite','positive','numel',2});
            validateattributes(fractions,{'numeric'},{'real','finite','nonnegative','numel',2});
            if sum(fractions)<=0, error('PSD:Fractions','Mixture fractions must have a positive sum.'); end
            fractions=fractions/sum(fractions);
            componentDensities=zeros(numel(diameter),2);
            for j=1:2
                componentDensities(:,j)=normalize_log_density(diameter,-log(diameter) ...
                    -0.5*(log(diameter/centers(j))/sigmas(j)).^2);
            end
            density=componentDensities*fractions;
            name='Bimodal';
    end
    distribution=struct('type',type,'name',name,'parameters',parameters, ...
        'diameter',diameter,'density',density,'cdf',cumtrapz(diameter,density), ...
        'componentDensities',componentDensities,'componentFractions',fractions, ...
        'effectiveParameters',effectiveParameters);
end

function value=positive_parameter(parameters,name)
    if ~isfield(parameters,name), error('PSD:MissingParameter','Missing parameter %s.',name); end
    value=parameters.(name);
    validateattributes(value,{'numeric'},{'scalar','real','finite','positive'});
end

function density=normalize_log_density(diameter,logDensity)
    offset=max(logDensity);
    if ~isfinite(offset), error('PSD:Underresolved','Distribution cannot be represented on this grid.'); end
    density=exp(logDensity-offset);
    density=density/trapz(diameter,density);
end
