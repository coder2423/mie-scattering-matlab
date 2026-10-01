function simulation = main_mie_scattering(overrides)
%MAIN_MIE_SCATTERING Run PSD generation, scattering comparison, plots and saving.
% Click Run, or call simulation=main_mie_scattering(). Default: preview, no files.
% Edit programs/default_scattering_config.m or pass a scalar overrides struct.
% Overrides replace top-level fields. For nested edits, modify a full default config.
    folder=fileparts(mfilename('fullpath')); addpath(fullfile(folder,'programs'));

    % -------------------- Load editable simulation parameters --------------------
    config=default_scattering_config();
    if nargin>=1 && ~isempty(overrides)
        if ~isstruct(overrides) || ~isscalar(overrides)
            error('Simulation:Options','Overrides must be a scalar struct.');
        end
        names=fieldnames(overrides);
        for j=1:numel(names)
            if ~isfield(config,names{j}), error('Simulation:UnknownOption','Unknown option %s.',names{j}); end
            config.(names{j})=overrides.(names{j});
        end
    end
    config=validate_scattering_config(config);
    simulation=run_scattering_simulation(config);
    nTypes=numel(simulation.distributions);

    % -------------------- PSD, forward-peak and scattering-tail plots --------------------
    if config.makePlots
        for w=1:numel(config.wavelengths_nm)
            fig=figure('Color','w','Name',sprintf('Mie and Fraunhofer: %g nm',config.wavelengths_nm(w)), ...
                'Visible',config.figureVisible,'Position',[80,50,1450,280*nTypes],'ToolBar','none');
            layout=tiledlayout(fig,nTypes,3,'Padding','compact','TileSpacing','compact');
            title(layout,sprintf('%g nm | particle index %.4g%+.3gi | host index %.4g', ...
                config.wavelengths_nm(w),real(config.nParticle),imag(config.nParticle),config.nMedium));
            targetAxes=gobjects(nTypes,3);
            for j=1:nTypes
                for k=1:3, targetAxes(j,k)=nexttile(layout); end
            end
            plot_scattering_comparison(simulation,w,targetAxes);
            if config.saveFigures
                if ~isfolder(config.outputFolder), mkdir(config.outputFolder); end
                name=sprintf('comparison_%.17gnm',config.wavelengths_nm(w));
                exportgraphics(fig,fullfile(config.outputFolder,[name,'.png']),'Resolution',160);
                savefig(fig,fullfile(config.outputFolder,[name,'.fig']));
            end

            % -------------------- Single-size, narrow-PSD and broad-PSD diagnostics --------------------
            if config.showPeakDiagnostics, plot_peak_diagnostics(simulation,w); end
        end
    end

    % -------------------- Save data only when explicitly enabled --------------------
    if config.saveData, save_scattering_data(simulation); end
    if config.saveFigures, fprintf('Figures saved to: %s\n',config.outputFolder); end
end
