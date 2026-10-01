function plot_peak_diagnostics(s,w)
%PLOT_PEAK_DIAGNOSTICS Show main peaks and side-lobe zooms of the same data.
    c=s.config;
    fig=figure('Color','w','Name',sprintf('Main-peak and side-lobe diagnostics: %g nm',c.wavelengths_nm(w)), ...
        'Visible',c.figureVisible,'Position',[100,50,1450,850],'ToolBar','none');
    layout=tiledlayout(fig,3,3,'Padding','compact','TileSpacing','compact');
    title(layout,sprintf('%g nm | Single size vs. narrow/broad PSDs: averaging of side lobes',c.wavelengths_nm(w)), ...
        'FontName','Helvetica');
    x=pi*c.nMedium*c.referenceDiameter_um*1e-6/(c.wavelengths_nm(w)*1e-9);
    if strcmp(c.fraunhoferVariant,'sphere'), firstMinimum=asin(3.8317059702075125/x)*180/pi;
    else, firstMinimum=3.8317059702075125/x*180/pi; end
    for j=1:3
        d=s.peakDiagnostics{j,w};
        % -------------------- Column 1: single-size or original narrow/broad PSD --------------------
        ax=nexttile(layout);
        if strcmp(d.sampling,'bin')
            stem(ax,d.diameter*1e6,1,'filled','LineWidth',1.5);
            ylabel(ax,'Single-size number fraction'); ylim(ax,[0,1.1]);
        else
            plot(ax,d.diameter*1e6,d.psd/1e6,'LineWidth',1.5);
            if strcmp(c.distributionBasis,'volume'), basisLabel='Volume'; else, basisLabel='Number'; end
            ylabel(ax,[basisLabel,' density (\mum^{-1})']);
        end
        xlim(ax,c.diameterRange_um); xlabel(ax,'Diameter (\mum)'); title(ax,d.name); grid(ax,'on');
        % -------------------- Column 2: independently peak-normalized main peaks --------------------
        axMain=nexttile(layout); hold(axMain,'on');
        plot(axMain,s.thetaDeg,d.mieNormalized,'LineWidth',1.5);
        plot(axMain,s.thetaDeg,d.fraunhoferNormalized,'--','LineWidth',1.5);
        xlim(axMain,[0,c.forwardMax_deg]); ylim(axMain,[0,1.05]);
        title(axMain,'Main peak: own-peak normalization'); xlabel(axMain,'Scattering angle (deg)');
        ylabel(axMain,'Normalized intensity'); grid(axMain,'on'); legend(axMain,{'Mie','Fraunhofer'});
        % -------------------- Column 3: side-lobe zoom of the same normalized curves --------------------
        axSide=nexttile(layout); hold(axSide,'on');
        plot(axSide,s.thetaDeg,d.mieNormalized,'LineWidth',1.5);
        plot(axSide,s.thetaDeg,d.fraunhoferNormalized,'--','LineWidth',1.5);
        xlim(axSide,[0,c.forwardMax_deg]); ylim(axSide,[0,c.sideLobeYMax]);
        if isreal(firstMinimum) && firstMinimum<c.forwardMax_deg
            xline(axSide,firstMinimum,':','First dark ring','HandleVisibility','off', ...
                'LabelVerticalAlignment','middle','FontName','Helvetica');
        end
        title(axSide,sprintf('Side lobes: y=0-%.3g; main-peak top clipped',c.sideLobeYMax));
        xlabel(axSide,'Scattering angle (deg)'); ylabel(axSide,'Normalized intensity'); grid(axSide,'on');
        axesList=[ax,axMain,axSide]; set(axesList,'FontName','Helvetica','FontSize',10,'XScale','linear');
        for a=axesList, a.Toolbar.Visible='off'; end
    end
    if c.saveFigures
        if ~isfolder(c.outputFolder), mkdir(c.outputFolder); end
        exportgraphics(fig,fullfile(c.outputFolder,sprintf('peak_diagnostics_%.17gnm.png',c.wavelengths_nm(w))),'Resolution',160);
        savefig(fig,fullfile(c.outputFolder,sprintf('peak_diagnostics_%.17gnm.fig',c.wavelengths_nm(w))));
    end
end
