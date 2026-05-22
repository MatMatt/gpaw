from gpaw.mixer.base import BaseMixer
from gpaw.mixer.pulay import PulayMixer
from gpaw.mixer.msr1 import MSR1Mixer
from gpaw.mixer.metric import BaseMetric, FFTMetric
# Imports for backwards compatibility:
from gpaw.old.mixer import _definemixerfunc

Mixer = _definemixerfunc('separate', 'pulay')
MixerSum = _definemixerfunc('sum', 'pulay')
MixerSum2 = _definemixerfunc('sum2', 'pulay')
MixerDif = _definemixerfunc('difference', 'pulay')
MixerFull = _definemixerfunc('fullspin', 'pulay')
FFTMixer = _definemixerfunc('separate', 'fft')
FFTMixerSum = _definemixerfunc('sum', 'fft')
FFTMixerSum2 = _definemixerfunc('sum2', 'fft')
FFTMixerDif = _definemixerfunc('difference', 'fft')
FFTMixerFull = _definemixerfunc('fullspin', 'fft')
BroydenMixer = _definemixerfunc('separate', 'broyden')
BroydenMixerSum = _definemixerfunc('sum', 'broyden')
BroydenMixerSum2 = _definemixerfunc('sum2', 'broyden')
BroydenMixerDif = _definemixerfunc('difference', 'broyden')


mixer_names = {
    'no-mixing': BaseMixer,
    'pulay': PulayMixer,
    'msr1': MSR1Mixer
}


def get_mixer_from_params(params: dict):
    # Ensure we don't touch the original dict:
    params = params.copy()

    # We should change how mixer metric is specified,
    # if we want to have more metric choices in the future.
    weight = params.pop('weight', 80)
    sigma = params.pop('sigma', 0.6)
    g_ss = params.pop('g_ss', None)
    if weight == 1:  # No metric
        metric = BaseMetric(g_ss=g_ss)
    else:
        metric = FFTMetric(g_ss=g_ss, weight=weight, sigma=sigma)
    params.pop('method', None)
    name = params.pop('backend', 'msr1')
    mixer = mixer_names[name](**params)
    mixer.metric = metric

    return mixer
