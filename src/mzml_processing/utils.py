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


def get_spectrum_collision_energy(spectrum: MSSpectrum) -> float:
    """Retrieves the collision energy of a spectrum, assuming structure
    <spectrum> <precursorList> <precursor> <activation> <cvParam name="collision energy" value=XX>
    and only one precursor entry.
    """
    return spectrum.getPrecursors()[0].getMetaValue("collision energy")
