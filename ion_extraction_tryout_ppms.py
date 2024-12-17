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
parser.add_argument(
    "--higher_collision_energy",
    type=int,
    help="Higher collision energy value",
    default=50.0,
)

args = parser.parse_args()

config = Config.from_path(Path("tooling_configs/ptm_search_experiment.json"))

exp = MSExperiment()
print(args.mzml_path)
MzMLFile().load(args.mzml_path, exp)

extractor = ScanWindowExtractor(
    config.lower_collision_energy, args.higher_collision_energy
)

exp_higher_energy = extractor.extract_higher_energy_windows(exp)

snr_threshold = 3

for ppm_tolerance in [5, 6, 7, 8, 9, 10]:
    detector = DiagnosticIonDetector(
        config.known_diagnostic_ions_file,
        ppm_tolerance,
        config.diagnostic_ions_mass_tolerance_unit,
        snr_threshold,
        args.higher_collision_energy,
    )

    print(datetime.datetime.now())
    print("Loaded spectra, starting ion extraction", flush=True)

    result_df = detector.extract_diagnostic_ions_for_spectra(exp_higher_energy.getSpectra())

    result_filename = f"{args.mzml_path}_diagnostic_ions_ppm_tolerance_{ppm_tolerance}_snr_threshold_{snr_threshold}_unimod.csv"

    print(datetime.datetime.now())
    print("Extracted spectra, saving to file", flush=True)
    result_df.to_csv(result_filename, index=False)
