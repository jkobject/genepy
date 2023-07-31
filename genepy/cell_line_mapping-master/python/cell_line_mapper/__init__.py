import requests
from collections import defaultdict, Counter
import pandas as pd
from io import StringIO

csv_url = 'https://intranet.broadinstitute.org/~datasci/cell_lines/name_mapping.csv'
__version__ = "0.1"

def read_file_to_dict(key_name, value_name):
    return_dict = defaultdict(str)

    response = requests.get(csv_url)
    assert response.status_code == 200, "Could not fetch mapping from {}, got, status_code={} reason={}".format(csv_url, response.status_code, response.reason)
    mapping_request = StringIO(response.text)
    df = pd.read_csv(mapping_request)

    for row in df.to_dict(orient='records'):
        if return_dict[row[key_name]] == "":
            return_dict[row[key_name]] = row[value_name]
        else:
            if row[value_name] != return_dict[row[key_name]]:      
                return_dict[row[key_name]] += str(", " + row[value_name])

    return return_dict


def check_unique(input_list, result_list):
    input_len = len(input_list)
    # make sure one key doesn't point to multiple values
    is_unique = True
    for mapping in result_list:
        if mapping is not None and ',' in mapping:
            is_unique = False
            break

    # make sure two keys don't point to the same value
    if is_unique:
        num_unique_values = len(Counter(tuple(mapping) for mapping in result_list if mapping is not None)) + int(None in result_list)

        is_unique = (input_len == num_unique_values)

    if not is_unique:
        raise RuntimeError('Mappings are not unique')


def name_mapper(input_type, input_names, output_type, ignore_problems=False, check_unique_mapping=True):
    if len(input_names) == 0:
        raise RuntimeError("Please input a non-empty list")

    mapping_dict = read_file_to_dict(input_type, output_type)
    output_names = []

    for i_name in input_names:
        o_name = mapping_dict.get(i_name)
        if not ignore_problems and o_name is None:
            raise KeyError(output_type + " " + "could not be found for the following " + input_type + ":  " + i_name)
        output_names.append(o_name)

    if check_unique_mapping:
        check_unique(input_names, output_names)

    assert len(output_names) == len(input_names)
    return output_names


def arxspan_to_ccle(arxspan_ids, ignore_problems=False, check_unique_mapping=True):
    return name_mapper('broad_id', arxspan_ids, 'canonical_ccle_name', ignore_problems, check_unique_mapping)


def ccle_to_arxspan(ccle_names, ignore_problems=False, check_unique_mapping=True):
    return name_mapper('ccle_name', ccle_names, 'broad_id', ignore_problems, check_unique_mapping)


def ccle_to_latest(ccle_names, ignore_problems=False, check_unique_mapping=True):
    return name_mapper('ccle_name', ccle_names, 'canonical_ccle_name', ignore_problems, check_unique_mapping)

# alias for the old name of this function
latest_ccle_names=ccle_to_latest
