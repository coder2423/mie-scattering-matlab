function report=test_scattering_output_controls()
%TEST_SCATTERING_OUTPUT_CONTROLS Verify preview and four independent save combinations.
% Temporary directories and a coarse grid test I/O behavior, not numerical accuracy.
    root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
    outputRoot=tempname;
    assert(startsWith(outputRoot,tempdir),'Test output must remain inside the system temporary directory.');
    mkdir(outputRoot);
    cleanup=onCleanup(@()rmdir(outputRoot,'s'));
    combinations=logical([0,0;1,0;0,1;1,1]);
    for j=1:size(combinations,1)
        destination=fullfile(outputRoot,sprintf('case_%d',j));
        before=findall(groot,'Type','figure');
        options=struct('wavelengths_nm',450,'distributionTypes',{{'normal'}}, ...
            'diameterPoints',31,'anglePoints',121,'diagnosticStd_um',[20,30], ...
            'figureVisible','off','outputFolder',destination);
        % Omit save flags in the first case to verify the preview-only default.
        if j>1
            options.saveData=combinations(j,1); options.saveFigures=combinations(j,2);
        end
        s=main_mie_scattering(options);
        figures=setdiff(findall(groot,'Type','figure'),before);
        assert(numel(figures)==2,'Expected comparison and diagnostic figures.');
        assert(s.config.saveData==combinations(j,1) && s.config.saveFigures==combinations(j,2));
        assert(isfile(fullfile(destination,'simulation_results.mat'))==combinations(j,1));
        assert(isfile(fullfile(destination,'comparison_metrics.csv'))==combinations(j,1));
        assert(isfile(fullfile(destination,'peak_diagnostic_narrow_450nm.csv'))==combinations(j,1));
        assert(isfile(fullfile(destination,'comparison_450nm.png'))==combinations(j,2));
        assert(isfile(fullfile(destination,'comparison_450nm.fig'))==combinations(j,2));
        assert(isfile(fullfile(destination,'peak_diagnostics_450nm.png'))==combinations(j,2));
        assert(isfile(fullfile(destination,'peak_diagnostics_450nm.fig'))==combinations(j,2));
        if j==1, assert(~isfolder(destination),'Preview must not create an output folder.'); end
        if combinations(j,1)
            loaded=load(fullfile(destination,'simulation_results.mat'),'simulation');
            assert(isequal(loaded.simulation.mieResults{1}.relativeIntensity,s.mieResults{1}.relativeIntensity));
        end
        close(figures);
    end
    report=struct('combinationCount',4,'outputRoot',outputRoot);
    disp(report); fprintf('SCATTERING_OUTPUT_CONTROLS_OK\n');
end
