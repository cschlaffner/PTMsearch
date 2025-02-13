from typing import Dict, FrozenSet, List, Union

from src.diagnostic_ions.utils import modification_unimod_format_to_dia_nn_varmod_format


def make_dia_nn_var_mod_commands(
    mods: List[Union[str, FrozenSet[str]]], additional_mods: List[str]
) -> List[str]:
    # add all modifications from the combination (or single mod)
    if mods == "unmodified":
        var_mod_commands = []
    else:
        mods_list = [mods] if isinstance(mods, str) else sorted(mods)
        var_mod_commands = [
            f"--var-mod {modification_unimod_format_to_dia_nn_varmod_format(mod)}"
            for mod in mods_list
        ]

    # add all mods that should be searched additionally
    var_mod_commands += [
        f"--var-mod {modification_unimod_format_to_dia_nn_varmod_format(mod)}"
        for mod in sorted(additional_mods)
    ]

    return var_mod_commands


def make_dia_nn_additional_param_commands(
    params: Dict[str, Union[int, str]]
) -> List[str]:
    return [f"--{param} {value}" for param, value in params.items()]
