function report = check_scattering_grid_convergence(mode,outputFolder)
%CHECK_SCATTERING_GRID_CONVERGENCE Refine diameter or angle sampling of main.
% Run main_mie_scattering with saveData=true first. Curve errors and cone-integral changes are
% numerical sampling differences, not rigorous error bounds of either model.
    folder=fileparts(fileparts(mfilename('fullpath')));
    if nargin<1 || isempty(mode), mode='diameter'; end
    if nargin<2 || isempty(outputFolder), outputFolder=fullfile(folder,'results','model_comparison'); end
    mode=validatestring(mode,{'diameter','angle'});
    file=fullfile(outputFolder,'simulation_results.mat');
    if ~isfile(file), error('Mie:MissingDemo','Run main_mie_scattering first.'); end
    loaded=load(file,'simulation'); s=loaded.simulation; c=s.config;
    d=s.diameter; theta=s.theta;
    if strcmp(mode,'diameter'), d=linspace(d(1),d(end),2*numel(d)-1).'; end
    if strcmp(mode,'angle'), theta=linspace(theta(1),theta(end),2*numel(theta)-1).'; end
    options=s.mieResults{1,1}.options;
    % Linear polarization is pi-periodic in azimuth. Resample its unwrapped
    % direction when the angle grid changes; a scalar direction is unchanged.
    if strcmp(mode,'angle') && ~isscalar(options.Azimuth)
        options.Azimuth=interp1(s.theta,unwrap(2*options.Azimuth(:))/2,theta,'linear');
    end
    rows=cell(numel(s.distributions)*numel(c.wavelengths_nm),1);
    for w=1:numel(c.wavelengths_nm)
        cached=[]; cachedF=[];
        for j=1:numel(s.distributions)
            g=generate_particle_distribution(d,s.distributions(j).type,s.distributions(j).parameters);
            m=mie_distribution_forward(theta,d,g.density,c.wavelengths_nm(w)*1e-9, ...
                c.nParticle,c.nMedium,options,cached);
            if j==1, cached=m; end
            f=fraunhofer_distribution_forward(theta,d,g.density,c.wavelengths_nm(w)*1e-9, ...
                c.nMedium,options,c.fraunhoferVariant,cachedF);
            if j==1, cachedF=struct('kernel',f.kernel,'cacheKey',f.cacheKey,'info',f.info); end
            m0=interp1(s.theta,s.mieResults{j,w}.differentialCrossSection,theta);
            f0=interp1(s.theta,s.fraunhoferResults{j,w}.differentialCrossSection,theta);
            mPeak=max(m.differentialCrossSection); fPeak=max(f.differentialCrossSection);
            if mPeak<=0 || fPeak<=0, error('Comparison:ZeroCurve','Cannot compare a zero curve.'); end
            diffM=abs(m0-m.differentialCrossSection); diffF=abs(f0-f.differentialCrossSection);
            metrics=compare_scattering_models(theta,m.differentialCrossSection,f.differentialCrossSection,c.forwardMax_deg*pi/180);
            previous=s.metrics((w-1)*numel(s.distributions)+j,:);
            rows{(w-1)*numel(s.distributions)+j}=struct('distribution',string(g.type), ...
                'wavelength_nm',c.wavelengths_nm(w),'miePeakScaledError',max(diffM)/mPeak, ...
                'mieMaxRelativeAboveFloor',max(diffM./max(m.differentialCrossSection,1e-8*mPeak)), ...
                'fraunhoferPeakScaledError',max(diffF)/fPeak, ...
                'mieConeIntegralRelativeChange',abs(metrics.mieFullConeCrossSection/previous.mieFullConeCrossSection-1), ...
                'modelAngularL2PercentRefined',metrics.angularL2Percent, ...
                'modelAngularL2ChangePercentagePoints',metrics.angularL2Percent-previous.angularL2Percent);
        end
        fprintf('GRID_%s %g nm completed\n',upper(mode),c.wavelengths_nm(w));
    end
    report=struct2table(vertcat(rows{:})); disp(report);
    writetable(report,fullfile(outputFolder,sprintf('grid_convergence_%s.csv',mode)));
    save(fullfile(outputFolder,sprintf('grid_convergence_%s.mat',mode)),'report');
end
