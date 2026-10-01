function save_scattering_data(simulation)
%SAVE_SCATTERING_DATA Save existing results as CSV and MAT without recomputing.
% Call save_scattering_data(simulation) explicitly after inspecting the figures.
    config=simulation.config;
    metrics=simulation.metrics; diameter=simulation.diameter;
    thetaDeg=simulation.thetaDeg; distributions=simulation.distributions;
    mieResults=simulation.mieResults; fraunhoferResults=simulation.fraunhoferResults;
    peakDiagnostics=simulation.peakDiagnostics;
    types=cellstr(string(config.distributionTypes));
    nTypes=numel(types); nWavelengths=numel(config.wavelengths_nm);
    if ~isfolder(config.outputFolder), mkdir(config.outputFolder); end
    writetable(metrics,fullfile(config.outputFolder,'comparison_metrics.csv'),'Encoding','UTF-8');
    for j=1:nTypes
        writetable(table(diameter*1e6,distributions(j).density/1e6,distributions(j).cdf, ...
            'VariableNames',{'diameter_um','density_per_um','cumulative_fraction'}), ...
            fullfile(config.outputFolder,sprintf('psd_%s.csv',types{j})));
        for w=1:nWavelengths
            mie=mieResults{j,w}; fraunhofer=fraunhoferResults{j,w};
            scale=max(mie.differentialCrossSection); if scale==0, scale=1; end
            exported=table(thetaDeg,mie.differentialCrossSection,fraunhofer.differentialCrossSection, ...
                mie.differentialCrossSection/scale,fraunhofer.differentialCrossSection/scale, ...
                mie.relativeIntensity,fraunhofer.relativeIntensity, ...
                'VariableNames',{'theta_deg','mie_dcs_m2_per_sr','fraunhofer_dcs_m2_per_sr', ...
                'mie_common_peak','fraunhofer_common_peak','mie_own_peak','fraunhofer_own_peak'});
            if ~isempty(config.totalNumber)
                exported.mie_irradiance_W_per_m2=mie.irradiance;
                exported.fraunhofer_irradiance_W_per_m2=fraunhofer.irradiance;
            end
            writetable(exported,fullfile(config.outputFolder,sprintf('angular_%s_%.17gnm.csv',types{j},config.wavelengths_nm(w))));
        end
    end
    if config.showPeakDiagnostics
        for w=1:nWavelengths
            for j=1:3
                d=peakDiagnostics{j,w};
                writetable(table(thetaDeg,d.mieNormalized,d.fraunhoferNormalized, ...
                    'VariableNames',{'theta_deg','mie_normalized_intensity','fraunhofer_normalized_intensity'}), ...
                    fullfile(config.outputFolder,sprintf('peak_diagnostic_%s_%.17gnm.csv',d.type,config.wavelengths_nm(w))));
            end
        end
    end
    save(fullfile(config.outputFolder,'simulation_results.mat'),'simulation','-v7.3');
    fprintf('Data saved to: %s\n',config.outputFolder);
end
