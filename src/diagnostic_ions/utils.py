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
