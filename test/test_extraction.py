from test.spectra import (
    COLLISION_ENERGY_HIGHER,
    COLLISION_ENERGY_LOWER,
    spectrum_ms1,
    spectrum_ms2_higher_energy,
    spectrum_ms2_lower_energy,
)
from typing import List

import pytest
from _pytest.fixtures import SubRequest
from pyopenms import MSExperiment, MSSpectrum, Precursor

from src.mzml_processing.extraction import ScanWindowExtractor


def spectra_ms1() -> List[MSSpectrum]:
    return [spectrum_ms1() for _ in range(2)]


def spectra_lower_energy() -> List[MSSpectrum]:
    return [spectrum_ms2_lower_energy() for _ in range(4)]


def spectra_higher_energy() -> List[MSSpectrum]:
    return [spectrum_ms2_higher_energy() for _ in range(3)]


@pytest.fixture
def experiment_higher_energy() -> MSExperiment:
    exp = MSExperiment()
    exp.setSpectra(spectra_higher_energy())
    return exp


@pytest.fixture
def experiment_empty() -> MSExperiment:
    return MSExperiment()


@pytest.fixture
def experiment_higher_lower_energy_ms1() -> MSExperiment:
    exp = MSExperiment()
    exp.setSpectra(spectra_higher_energy() + spectra_lower_energy() + spectra_ms1())
    return exp


@pytest.fixture
def experiment_lower_energy_ms1() -> MSExperiment:
    exp = MSExperiment()
    exp.setSpectra(spectra_lower_energy() + spectra_ms1())
    return exp


@pytest.fixture
def experiment_higher_energy_ms1() -> MSExperiment:
    exp = MSExperiment()
    exp.setSpectra(spectra_higher_energy() + spectra_ms1())
    return exp


@pytest.fixture
def experiment_ms1() -> MSExperiment:
    exp = MSExperiment()
    exp.setSpectra(spectra_ms1())
    return exp


@pytest.mark.parametrize(
    "exp_input_fixture, exp_output_fixture",
    [
        ("experiment_empty", "experiment_empty"),
        ("experiment_higher_energy", "experiment_empty"),
        ("experiment_ms1", "experiment_ms1"),
        ("experiment_higher_energy_ms1", "experiment_ms1"),
        ("experiment_lower_energy_ms1", "experiment_lower_energy_ms1"),
        ("experiment_higher_lower_energy_ms1", "experiment_lower_energy_ms1"),
    ],
)
def test_extract_lower_energy(
    exp_input_fixture: str,
    exp_output_fixture: str,
    request: SubRequest,
) -> None:
    exp_input = request.getfixturevalue(exp_input_fixture)
    exp_output = request.getfixturevalue(exp_output_fixture)

    extractor = ScanWindowExtractor(
        COLLISION_ENERGY_LOWER,
        COLLISION_ENERGY_HIGHER,
    )
    result = extractor.extract_ms1_and_lower_energy_windows(exp_input)
    assert result == exp_output
