function [numberFractions,weights,numberScale] = particle_distribution_weights(diameter,psd,basis,sampling)
%PARTICLE_DISTRIBUTION_WEIGHTS Shared PSD quadrature for both optical models.
% Density is per meter of diameter; bin inputs are already integrated fractions.
    basis=validatestring(basis,{'number','volume'});
    sampling=validatestring(sampling,{'density','bin'});
    validateattributes(diameter,{'numeric'},{'vector','nonempty','real','finite','positive'});
    validateattributes(psd,{'numeric'},{'vector','nonempty','real','finite','nonnegative'});
    diameter=diameter(:); psd=psd(:);
    if numel(diameter)~=numel(psd) || any(diff(diameter)<=0)
        error('Mie:Grid','diameter must increase strictly and match psd length.');
    end
    weights=ones(size(diameter));
    if strcmp(sampling,'density')
        if numel(diameter)<2
            error('Mie:DensityGrid','Density needs >=2 nodes; use bin for monodisperse input.');
        end
        spacing=diff(diameter);
        weights=[spacing(1);spacing(1:end-1)+spacing(2:end);spacing(end)]/2;
    end
    if strcmp(basis,'volume'), weights=weights./(pi*diameter.^3/6); end
    raw=weights.*psd;
    numberScale=sum(raw);
    if ~isfinite(numberScale) || numberScale<=0
        error('Mie:ZeroPSD','PSD must have a finite positive integral/total.');
    end
    numberFractions=raw/numberScale;
end
