"""Generation of MSSpectrum instances for testing."""

from typing import Optional, Tuple

import numpy as np
from pyopenms import MSSpectrum, Precursor

COLLISION_ENERGY_HIGHER = 50
COLLISION_ENERGY_LOWER = 30


def spectrum(
    ms_level: int,
    native_id: Optional[str] = None,
    collision_energy: Optional[int] = None,
    peaks_mz_intensities: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> MSSpectrum:
    spec = MSSpectrum()
    spec.setMSLevel(ms_level)

    if native_id is not None:
        spec.setNativeID(native_id)

    if collision_energy is not None:
        precursor = Precursor()
        precursor.setMetaValue("collision energy", collision_energy)
        spec.setPrecursors([precursor])

    if peaks_mz_intensities is not None:
        spec.set_peaks(peaks_mz_intensities)

    return spec


def spectrum_ms1(native_id: Optional[str] = None) -> MSSpectrum:
    """Generates an MS1 spectrum. Collision energy does not apply here.
    M/z and intensities of MS1 spectra are not relevant in this software
    at the current and planned state."""
    return spectrum(1, native_id=native_id)


def spectrum_ms2_lower_energy(
    native_id: Optional[str] = None,
) -> MSSpectrum:
    """Generates an MS2 spectrum with lower collision energy. M/z and intensities
    of lower-energy spectra are not relevant in this software
    at the current and planned state."""
    return spectrum(
        2,
        native_id=native_id,
        collision_energy=COLLISION_ENERGY_LOWER,
    )


def spectrum_ms2_higher_energy(
    native_id: Optional[str] = None,
    peaks_mz_intensities: Optional[Tuple[np.ndarray, np.ndarray]] = None,
) -> MSSpectrum:
    """Generates an MS2 spectrum with higher collision energy."""
    return spectrum(
        2,
        native_id=native_id,
        collision_energy=COLLISION_ENERGY_HIGHER,
        peaks_mz_intensities=peaks_mz_intensities,
    )
