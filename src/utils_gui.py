from typing import Any, Dict, FrozenSet, List, Union


def get_mod_unimod(mod_input: str) -> str:
    amino_acid = mod_input[0]
    unimod_id = mod_input[1:]
    return f"{amino_acid}(UniMod:{unimod_id})"


def get_mod_list(mod_list_input: str) -> List[str]:
    if mod_list_input == "":
        return []
    return [get_mod_unimod(mod) for mod in mod_list_input.split(";")]


def get_mod_combinations_list(mod_combinations_list_input: str) -> List[FrozenSet[str]]:
    if mod_combinations_list_input == "":
        return []
    mod_combinations_list = mod_combinations_list_input.split(";")
    mod_combinations_sets_list = []
    for combination in mod_combinations_list:
        mod_list = [get_mod_unimod(mod) for mod in combination.split("|")]
        mod_combinations_sets_list.append(frozenset(mod_list))
    return mod_combinations_sets_list


def get_dictionary(dict_input: str) -> Dict[str, Any]:
    if dict_input == "":
        return {}

    def try_to_number(value):
        try:
            value_float = float(value)
            if value_float % 1 == 0:
                return int(value_float)
            return value_float
        except (ValueError, TypeError):
            return value

    dict_list = dict_input.split(";")
    key_value_list = [key_value.split(":") for key_value in dict_list]
    dict_config = {
        key_value[0]: try_to_number(key_value[1]) for key_value in key_value_list
    }
    return dict_config


def get_dictionary_libraries_by_mod(
    dict_input: str,
) -> Dict[Union[str, FrozenSet], str]:
    if dict_input == "":
        return {}

    libraries_dict = get_dictionary(dict_input)
    libraries_dict_corrected_keys = {}
    for mod_string, library_path in libraries_dict.items():
        mod_key = (
            frozenset([get_mod_unimod(mod) for mod in mod_string.split("|")])
            if "|" in mod_string
            else get_mod_unimod(mod_string)
        )
        libraries_dict_corrected_keys[mod_key] = library_path
    return libraries_dict_corrected_keys
