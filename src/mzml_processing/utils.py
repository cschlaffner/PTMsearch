from typing import Union, cast

from pyopenms import MSSpectrum, MzMLFile, PeakFileOptions


def get_diann_compatible_mzml_output_file() -> MzMLFile:
    """Creates an MzMLFile with options set based on
    https://github.com/vdemichev/DiaNN?tab=readme-ov-file#raw-data-formats
    (i.e., 32-bit encoding, no compression)"""
    options = PeakFileOptions()
    options.setMz32Bit(True)
    options.setIntensity32Bit(True)
    options.setCompression(False)

    output_file = MzMLFile()
    output_file.setOptions(options)
    return output_file


def check_ms1_spectrum(spectrum: MSSpectrum) -> bool:
    """Checks whether a spectrum is of MS level 1."""
    return spectrum.getMSLevel() == 1


def validate_ms1_spectrum(spectrum: MSSpectrum) -> bool:
    """Validates that spectrum is of MS level 1."""
    assert check_ms1_spectrum(spectrum), (
        f"Spectrum {spectrum.getNativeID()} should be an MS1 scan, "
        f"but has MSLevel {spectrum.getMSLevel()}."
    )


def check_ms2_spectrum(spectrum: MSSpectrum) -> bool:
    """Checks whether a spectrum is of MS level 2."""
    return spectrum.getMSLevel() == 2


def validate_ms2_spectrum(spectrum: MSSpectrum) -> bool:
    """Validates that spectrum is of MS level 2."""
    assert check_ms2_spectrum(spectrum), (
        f"Spectrum {spectrum.getNativeID()} should be an MS2 scan, "
        f"but has MSLevel {spectrum.getMSLevel()}"
    )


def get_spectrum_collision_energy(spectrum: MSSpectrum) -> Union[float, int]:
    """Retrieves the collision energy of a spectrum, assuming structure
    <spectrum> <precursorList> <precursor> <activation> <cvParam name="collision energy" value=XX>
    and only one precursor entry.
    """
    validate_ms2_spectrum(spectrum)

    precursors = spectrum.getPrecursors()
    assert len(precursors) > 0, (
        f"Spectrum {spectrum.getNativeID()} has no precursors to extract "
        "the collision energy from, check your MzML structure."
    )
    collision_energy = precursors[0].getMetaValue("collision energy")
    assert isinstance(
        collision_energy, (float, int)
    ), "Collision energy value could not be extracted, check your MzML structure."

    # Type casting to please the typecheck
    return cast(float, collision_energy)


def check_collision_energy(
    spectrum: MSSpectrum, collision_energy: Union[float, int]
) -> bool:
    """Checks whether a spectrum has a certain collision energy."""
    return get_spectrum_collision_energy(spectrum) == collision_energy


def check_collision_energy_ms2_spectrum(
    spectrum: MSSpectrum, collision_energy: Union[float, int]
) -> bool:
    """Checks whether a spectrum is of MS level 2 and with a certain collision energy."""
    return check_ms2_spectrum(spectrum) and check_collision_energy(
        spectrum, collision_energy
    )


def validate_collision_energy_ms2_spectrum(
    spectrum: MSSpectrum, collision_energy: Union[float, int]
) -> None:
    """Validates that a spectrum is of MS level 2 and with a certain collision energy."""
    validate_ms2_spectrum(spectrum)

    assert check_collision_energy(spectrum, collision_energy), (
        f"Spectrum {spectrum.getNativeID()} should have collision energy {collision_energy} "
        f"but has {get_spectrum_collision_energy(spectrum)}."
    )
