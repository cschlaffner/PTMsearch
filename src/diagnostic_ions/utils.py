import re

from pyopenms import ModificationsDB, ResidueDB

modifications_db = ModificationsDB()
residue_db = ResidueDB()


def get_modification_unimod_format(amino_acid: str, mod_name: str) -> str:
    """Convert a modification given as amino acid (full name) and modification
    (PSI-MS Name (if available) or Interim Name) to the UniMod format that is used
    in spectral libraries predicted by DIA-NN.
    E.g. Tyrosine,Phospho -> Y(UniMod:21)
    """
    unimod_accession = modifications_db.getModification(mod_name).getUniModAccession()
    amino_acid_letter = residue_db.getResidue(amino_acid).getOneLetterCode()
    return f"{amino_acid_letter}({unimod_accession})"


def modification_unimod_format_to_dia_nn_varmod_format(modification: str) -> str:
    """Convert a modification in UniMod format to the format used in the --var-mod
    command of DIA-NN, including monoisotopic mass.
    E.g. Y(UniMod:21) -> UniMod:21,79.966331,Y"""
    amino_acid = modification[0]
    mod_unimod = re.findall("UniMod:[0-9]+", modification)[0]
    mod_mass = modifications_db.getModification(mod_unimod).getDiffMonoMass()
    return f"{mod_unimod},{mod_mass},{amino_acid}"


def diff_mono_mass_to_unimod_format(
    amino_acid_letter: str, diff_mono_mass: float
) -> str:
    """Convert a modification given as amino acid (letter) and
    monoisotopic mass difference to the UniMod format that is used
    in spectral libraries predicted by DIA-NN.
    E.g. Y, 79.966331 -> Y(UniMod:21)
    """
    unimod_accession = modifications_db.getBestModificationByDiffMonoMass(
        diff_mono_mass, 0.0001, amino_acid_letter, 0
    ).getUniModAccession()
    return f"{amino_acid_letter}({unimod_accession})"
