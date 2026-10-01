function [yM,yF,label]=scattering_plot_values(mie,fraunhofer,quantity)
%SCATTERING_PLOT_VALUES Select display units without changing physical results.
    switch quantity
        case 'own_peak'
            yM=mie.relativeIntensity; yF=fraunhofer.relativeIntensity;
            label='Normalized intensity I / I_{max}';
        case 'common_peak'
            scale=max(mie.differentialCrossSection); if scale==0, scale=1; end
            yM=mie.differentialCrossSection/scale; yF=fraunhofer.differentialCrossSection/scale;
            label='Intensity / Mie peak';
        case 'dcs'
            yM=mie.differentialCrossSection; yF=fraunhofer.differentialCrossSection; label='Mean differential cross section (m^2/sr)';
        case 'irradiance'
            if isempty(mie.irradiance) || isempty(fraunhofer.irradiance)
                error('Simulation:Irradiance','Set TotalNumber before requesting irradiance.');
            end
            yM=mie.irradiance; yF=fraunhofer.irradiance; label='Scattered irradiance (W/m^2)';
        otherwise
            error('Simulation:PlotQuantity','Unknown plot quantity: %s.',quantity);
    end
end
