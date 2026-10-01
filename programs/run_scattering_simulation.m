function simulation = run_scattering_simulation(config)
%RUN_SCATTERING_SIMULATION Generate PSDs and compare models without plotting or saving.
    diameter=linspace(config.diameterRange_um(1),config.diameterRange_um(2),config.diameterPoints).'*1e-6;
    thetaDeg=linspace(config.angleRange_deg(1),config.angleRange_deg(2),config.anglePoints).';
    theta=thetaDeg*pi/180; wavelengthList=config.wavelengths_nm(:)*1e-9;
    types=cellstr(string(config.distributionTypes));
    for j=1:numel(types)
        if ~isfield(config.parameters,types{j}), error('Simulation:Distribution','Missing parameters for %s.',types{j}); end
        current=generate_particle_distribution(diameter,types{j},config.parameters.(types{j}));
        if j==1, distributions=current; else, distributions(j)=current; end
    end
    options=struct('DistributionBasis',config.distributionBasis,'DistributionSampling','density', ...
        'Polarization',config.polarization,'Azimuth',config.azimuth_rad, ...
        'IncidentIrradiance',config.incidentIrradiance,'Distance',config.distance_m,'TotalNumber',config.totalNumber);
    nTypes=numel(types); nWavelengths=numel(wavelengthList);
    mieResults=cell(nTypes,nWavelengths); fraunhoferResults=mieResults;
    mieKernels=cell(nWavelengths,1); fraunhoferKernels=mieKernels;
    elapsedSeconds=zeros(nWavelengths,1); metricsRows=cell(nTypes*nWavelengths,1);
    peakDiagnostics=cell(3,nWavelengths);

    % Build each wavelength kernel once and reuse it for all distributions.
    for w=1:nWavelengths
        started=tic; cached=[]; cachedF=[];
        for j=1:nTypes
            mie=mie_distribution_forward(theta,diameter,distributions(j).density, ...
                wavelengthList(w),config.nParticle,config.nMedium,options,cached);
            fraunhofer=fraunhofer_distribution_forward(theta,diameter,distributions(j).density, ...
                wavelengthList(w),config.nMedium,options,config.fraunhoferVariant,cachedF);
            if j==1
                mieKernels{w}=struct('kernel',mie.kernel,'cacheKey',mie.cacheKey, ...
                    'efficiencies',mie.efficiencies,'crossSections',mie.crossSections);
                cached=mieKernels{w};
                cachedF=struct('kernel',fraunhofer.kernel,'cacheKey',fraunhofer.cacheKey,'info',fraunhofer.info);
                fraunhoferKernels{w}=fraunhofer.kernel;
            end
            comparison=compare_scattering_models(theta,mie.differentialCrossSection, ...
                fraunhofer.differentialCrossSection,config.forwardMax_deg*pi/180);
            row=struct('distribution',string(types{j}),'wavelength_nm',config.wavelengths_nm(w));
            fields=fieldnames(comparison);
            for k=1:numel(fields), row.(fields{k})=comparison.(fields{k}); end
            metricsRows{(w-1)*nTypes+j}=row;
            % Store large matrices once, outside the per-distribution results.
            mieResults{j,w}=rmfield(mie,{'kernel','distributionKernel'});
            fraunhoferResults{j,w}=rmfield(fraunhofer,{'kernel','distributionKernel'});
        end
        if config.showPeakDiagnostics
            for j=1:3
                if j==1
                    diagnosticD=config.referenceDiameter_um*1e-6;
                    diagnosticPsd=1; diagnosticSampling='bin';
                    diagnosticName='Single size'; diagnosticType='single'; diagnosticCache=[]; diagnosticCacheF=[];
                else
                    diagnosticD=diameter;
                    g=generate_particle_distribution(diameter,'normal', ...
                        struct('Mean',config.referenceDiameter_um*1e-6,'Std',config.diagnosticStd_um(j-1)*1e-6));
                    diagnosticPsd=g.density; diagnosticSampling='density'; diagnosticCache=mieKernels{w}; diagnosticCacheF=cachedF;
                    if j==2, diagnosticType='narrow'; else, diagnosticType='broad'; end
                    diagnosticName=sprintf('Normal: standard deviation %g um',config.diagnosticStd_um(j-1));
                end
                diagnosticOptions=options; diagnosticOptions.DistributionSampling=diagnosticSampling;
                m=mie_distribution_forward(theta,diagnosticD,diagnosticPsd,wavelengthList(w), ...
                    config.nParticle,config.nMedium,diagnosticOptions,diagnosticCache);
                f=fraunhofer_distribution_forward(theta,diagnosticD,diagnosticPsd,wavelengthList(w), ...
                    config.nMedium,diagnosticOptions,config.fraunhoferVariant,diagnosticCacheF);
                peakDiagnostics{j,w}=struct('name',diagnosticName,'type',diagnosticType, ...
                    'diameter',diagnosticD,'psd',diagnosticPsd,'sampling',diagnosticSampling, ...
                    'mieNormalized',m.relativeIntensity,'fraunhoferNormalized',f.relativeIntensity, ...
                    'mieDcs',m.differentialCrossSection,'fraunhoferDcs',f.differentialCrossSection);
            end
        end
        elapsedSeconds(w)=toc(started);
        fprintf('%g nm: %d distributions, one Mie kernel, %.2f s\n',config.wavelengths_nm(w),nTypes,elapsedSeconds(w));
    end
    metrics=struct2table(vertcat(metricsRows{:}));
    simulation=struct('config',config,'theta',theta,'thetaDeg',thetaDeg,'diameter',diameter, ...
        'distributions',distributions,'mieResults',{mieResults},'fraunhoferResults',{fraunhoferResults}, ...
        'mieKernels',{mieKernels},'fraunhoferKernels',{fraunhoferKernels},'metrics',metrics, ...
        'peakDiagnostics',{peakDiagnostics},'elapsedSeconds',elapsedSeconds);
    disp(metrics(:,{'distribution','wavelength_nm','angularL2Percent', ...
        'independentPeakShapeL2Percent','forwardIntensityDifferencePercent', ...
        'forwardConeDifferencePercent','fullConeDifferencePercent'}));

end
