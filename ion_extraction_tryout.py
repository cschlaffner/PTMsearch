import datetime
from argparse import ArgumentParser
from pathlib import Path

from pyopenms import MSExperiment, MSSpectrum, MzMLFile

from src.config.config import Config
from src.diagnostic_ions.detection import DiagnosticIonDetector
from src.mzml_processing.extraction import ScanWindowExtractor

parser = ArgumentParser()
parser.add_argument(
    "mzml_path",
    type=str,
    help="Path to mzML file with higher and lower collision energy windows",
)

args = parser.parse_args()

config = Config.from_path(Path("tooling_configs/ptm_search_experiment.json"))

exp = MSExperiment()
MzMLFile().load(args.mzml_path, exp)

extractor = ScanWindowExtractor(
    config.lower_collision_energy, config.higher_collision_energy
)

exp_higher_energy = extractor.extract_higher_energy_windows(exp)

ppm_tolerance = 5

for snr_threshold in [1, 2, 3, 5, 10, 20, 50, 100]:

    detector = DiagnosticIonDetector(
        config.known_diagnostic_ions_file,
        ppm_tolerance,
        config.diagnostic_ions_mass_tolerance_unit,
        snr_threshold,
        config.higher_collision_energy,
    )

    print(datetime.datetime.now())
    print("Loaded spectra, starting ion extraction", flush=True)

    result_df = detector.extract_diagnostic_ions_for_spectra(
        exp_higher_energy.getSpectra()
    )

    result_filename = f"{args.mzml_path}_diagnostic_ions_ppm_tolerance_{ppm_tolerance}_snr_threshold_{snr_threshold}.csv"

    print(datetime.datetime.now())
    print("Extracted spectra, saving to file", flush=True)
    result_df.to_csv(result_filename, index=False)
