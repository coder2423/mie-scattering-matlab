function report=test_grid_convergence()
%TEST_GRID_CONVERGENCE Verify saved-run refinement with angle-dependent azimuth.
    root=fileparts(fileparts(mfilename('fullpath')));
    addpath(root,fullfile(root,'programs'));
    destination=tempname;
    assert(startsWith(destination,tempdir),'Test output must remain inside the system temporary directory.');
    cleanup=onCleanup(@()rmdir(destination,'s'));
    c=default_scattering_config();
    c.wavelengths_nm=450; c.distributionTypes={'normal'};
    c.diameterPoints=11; c.anglePoints=21; c.angleRange_deg=[0,2];
    c.makePlots=false; c.showPeakDiagnostics=false;
    c.polarization='linear'; c.azimuth_rad=linspace(0,0.4,c.anglePoints).';
    c.saveData=true; c.outputFolder=destination;
    main_mie_scattering(c);
    angle=check_scattering_grid_convergence('angle',destination);
    diameter=check_scattering_grid_convergence('diameter',destination);
    assert(height(angle)==1 && height(diameter)==1);
    assert(all(isfinite(angle{:,3:end}),'all') && all(isfinite(diameter{:,3:end}),'all'));
    report=struct('angleDependentAzimuthRefinement',true,'refinementModes',2);
    disp(report); fprintf('GRID_CONVERGENCE_WORKFLOW_OK\n');
end
