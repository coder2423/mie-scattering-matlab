function report=test_release_readiness()
%TEST_RELEASE_READINESS Exercise cached diffraction, preflight and public edge cases.
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(root,fullfile(root,'programs'));
    c=default_scattering_config();
    assert(strcmp(c.outputFolder,fullfile(root,'results','model_comparison')));
    assert(~c.saveData && ~c.saveFigures);
    assert(~contains(which('default_scattering_config'),[filesep,'apps',filesep]));

    d=linspace(10,200,31).'*1e-6; theta=linspace(0,2,121).'*pi/180;
    first=generate_particle_distribution(d,'normal',struct('Mean',100e-6,'Std',20e-6));
    second=generate_particle_distribution(d,'rr',struct('Mode',100e-6,'Shape',3));
    options=struct('DistributionBasis','volume');
    f=fraunhofer_distribution_forward(theta,d,first.density,450e-9,1,options);
    reused=fraunhofer_distribution_forward(theta,d,second.density,450e-9,1,options,'sphere',f);
    direct=fraunhofer_distribution_forward(theta,d,second.density,450e-9,1,options);
    assert(reused.cacheHit && ~direct.cacheHit);
    assert(isequal(reused.differentialCrossSection,direct.differentialCrossSection));
    options.TotalNumber=100; options.IncidentIrradiance=2; options.Distance=4;
    scaled=fraunhofer_distribution_forward(theta,d,second.density,450e-9,1,options,'sphere',f);
    assert(isequal(scaled.irradiance,12.5*direct.differentialCrossSection));
    must_error(@()fraunhofer_distribution_forward(theta,d,second.density,632e-9,1,options,'sphere',f),'Fraunhofer:CacheMismatch');
    must_error(@()fraunhofer_distribution_forward(theta,d,second.density,450e-9,1,options,'paraxial',f),'Fraunhofer:CacheMismatch');
    must_error(@()fraunhofer_distribution_forward(theta,d,second.density,450e-9,1,options,'sphere',struct()),'Fraunhofer:CacheMismatch');
    m=mie_distribution_forward(theta,d,first.density,450e-9,1.591,1);
    must_error(@()mie_distribution_forward(theta,d,first.density,450e-9,1.591,1,struct(),struct('cacheKey',m.cacheKey)),'Mie:CacheMismatch');

    must_error(@()main_mie_scattering(struct('wavelengths_nm',[450,450])),'Simulation:Wavelengths');
    must_error(@()main_mie_scattering(struct('nParticle',1.5+0.01i)),'Mie:IndexConvention');
    must_error(@()main_mie_scattering(struct('nParticle',1)),'Comparison:ZeroMie');
    must_error(@()main_mie_scattering(struct('outputFolder','')),'Simulation:OutputFolder');
    must_error(@()main_mie_scattering(struct('azimuth_rad',[0,1])),'Mie:Azimuth');
    must_error(@()main_mie_scattering(struct('distributionTypes',{{'unknown'}})),'Simulation:Types');
    must_error(@()main_mie_scattering(struct('diameterRange_um',[10,1e6])),'Mie:SizeRange');
    must_error(@()main_mie_scattering(struct('makePlots',false,'saveFigures',true)),'Simulation:FigureSaving');
    % Close decimal wavelengths must remain distinct in exported filenames.
    s=main_mie_scattering(struct('wavelengths_nm',[450,450.00001], ...
        'diameterPoints',7,'anglePoints',11,'showPeakDiagnostics',false,'makePlots',false));
    destination=tempname;
    assert(startsWith(destination,tempdir),'Test output must remain inside the system temporary directory.');
    mkdir(destination); cleanup=onCleanup(@()rmdir(destination,'s'));
    s.config.outputFolder=destination; save_scattering_data(s);
    assert(numel(dir(fullfile(destination,'angular_normal_*nm.csv')))==2);
    assert(all(cellfun(@(r)r.cacheHit,s.fraunhoferResults(2:end,1))));
    report=struct('diffractionCacheExactAgreement',true,'preflightCases',8,'distinctDecimalExports',true);
    disp(report); fprintf('RELEASE_READINESS_OK\n');
end

function must_error(action,identifier)
    try
        action();
    catch exception
        assert(strcmp(exception.identifier,identifier),'Unexpected error identifier: %s.',exception.identifier);
        return;
    end
    error('Test:ExpectedError','Expected error %s was not raised.',identifier);
end
