from . import name_mapper, arxspan_to_ccle, ccle_to_arxspan, latest_ccle_names, csv_url
import cell_line_mapper
import pytest
import pandas as pd
import requests
from io import StringIO
from collections import defaultdict
import csv
##################################################################################
### TEST ARXSPAN_TO_CCLE
#just the first three broadid's in the doc

def test_real_data_fetch_check():
    return_dict = defaultdict(set)

    mapping_request = StringIO(requests.get(csv_url).text)
    df = pd.read_csv(mapping_request)
    assert list(df) == ['ccle_name', 'canonical_ccle_name', 'broad_id']


##################################################################################

@pytest.fixture()
def fake_mapping_csv(monkeypatch):
    def mock_read_file_to_dict(key_name, value_name):
        return_dict = defaultdict(str)

        with open('test-data.csv', mode='r') as csvfile:
            map_reader = csv.DictReader(csvfile)
            for rows in map_reader:
                if return_dict[rows[key_name]] == "":
                    return_dict[rows[key_name]] = rows[value_name]
                else:
                    return_dict[rows[key_name]]+=str(", "+rows[value_name])
        print(return_dict)
        return return_dict

    monkeypatch.setattr(cell_line_mapper, 'read_file_to_dict', mock_read_file_to_dict)


##################################################################################
### TEST ARXSPAN_TO_CCLE
# just the first three broadid's in the doc
def test_arxspan_to_ccle_first_three_rows(fake_mapping_csv):
    assert arxspan_to_ccle(["1", "2", "3"]) == ["A", "B", "C"]


#  1.  unmappable ID and ignore_problems == false
def test_arxspan_to_ccle_unmappable_ID_ignore_problems_false(fake_mapping_csv):
    with pytest.raises(KeyError) as excinfo:
        arxspan_to_ccle(["1", "madeupfakename"])
    assert "canonical_ccle_name could not be found for the following broad_id:  madeupfakename" in str(excinfo.value)


#  2.  unmappable ID and ignore_problems == true
def test_arxspan_to_ccle_unmappable_ID_ignore_problems_true(fake_mapping_csv):
    assert arxspan_to_ccle(["1", "madeupfakename"], True) == ["A", None]


#  2.5.  arxspan id has multiple ccle names
def test_arxspan_to_ccle_one_key_many_values(fake_mapping_csv):
    assert arxspan_to_ccle(["7"], True, False) == ["D, F"]


#  3.  unique elements in arxspan_ids do not map to unique ccle names and check_unique_mapping == false
def test_arxspan_to_ccle_nonunique_mapping_check_unique_mapping_false(fake_mapping_csv):
    assert arxspan_to_ccle(["4", "5"], True, False) == ["D", "D"]


#  4.  unique elements in arxspan_ids do not map to unique ccle names and check_unique_mapping == true
def test_arxspan_to_ccle_nonunique_mapping_check_unique_mapping_true(fake_mapping_csv):
    with pytest.raises(RuntimeError) as excinfo:
        arxspan_to_ccle(["4", "5"], True, True)
    assert 'Mappings are not unique' in str(excinfo.value)


##################################################################################
### TEST CCLE_ARXSPAN
#  0.  Just the first three rows
def test_ccle_to_arxspan_first_three_rows(fake_mapping_csv):
    assert ccle_to_arxspan(["a", "B", "c"]) == ["1", "2", "3"]


#  1.  unmappable ccle name and ignore_problems == false
def test_ccle_to_arxspan_unmappable_ID_ignore_problems_false(fake_mapping_csv):
    with pytest.raises(KeyError) as excinfo:
        ccle_to_arxspan(["a", "madeupfakename", "idk"], False)
        assert "broad_id could not be found for the following canonical_ccle_name:  madeupfakename" in str(excinfo.value)


#  2.  unmappable ccle name and ignore_problems == true
def test_ccle_to_arxspan_unmappable_ID_ignore_problems_true(fake_mapping_csv):
    assert ccle_to_arxspan(["a", "madeupfakename"], True) == ["1", None]


#  2.5.  ccle name has multiple arxspan_ids
def test_ccle_to_arxspan_one_key_many_values(fake_mapping_csv):
    assert ccle_to_arxspan(["d"], True, False) == ["4, 5, 7"]


#  3.  unique ccle names do not map to unique arxspan_ids and check_unique_mapping == false
def test_ccle_to_arxspan_nonunique_mapping_check_unique_mapping_false(fake_mapping_csv):
    assert ccle_to_arxspan(["e"], True, False) == ["6, 7"]


#  4.  unique ccle names do not map to unique arxspan_ids and check_unique_mapping == true
def test_ccle_to_arxspan_nonunique_mapping_check_unique_mapping_true(fake_mapping_csv):
    with pytest.raises(RuntimeError) as excinfo:
        a = ccle_to_arxspan(["e"], True, True)
    assert 'Mappings are not unique' in str(excinfo.value)


##################################################################################
### TEST LATEST_CCLE_NAMES

#  0.  Just a few samples that have different names
def test_latest_ccle_names_first_three_rows(fake_mapping_csv):
    assert latest_ccle_names(["a", "B", "c"]) == ["A", "B", "C"]


#  1.  unmappable ccle name and ignore_problems == false
def test_latest_ccle_names_unmappable_ignore_problems_false(fake_mapping_csv):
    with pytest.raises(KeyError) as excinfo:
        latest_ccle_names(["madeupfakename"], False)
    assert "canonical_ccle_name could not be found for the following ccle_name:  madeupfakename" in str(excinfo.value)


#  2.  unmappable ccle name and ignore_problems == true
def test_latest_ccle_names_unmappable_ignore_problems_true(fake_mapping_csv):
    assert latest_ccle_names(["madeupfakename", "a"], True) == [None, "A"]


#    2.5  old ccle name maps to multiple latest names
def test_latest_ccle_names_one_key_many_values(fake_mapping_csv):
    assert latest_ccle_names(["e"], True, False) == ["E, F"]

def test_latest_ccle_names_one_key_many_values_check_unique_mapping_true(fake_mapping_csv):
    with pytest.raises(RuntimeError) as excinfo:
        latest_ccle_names(["e"], True, True)
    assert 'Mappings are not unique' in str(excinfo.value)

#  4.  unique ccle names do not map to unique latest names and check_unique_mapping == true
def test_latest_ccle_names_nonunique_mapping_check_unique_mapping_true(fake_mapping_csv):
    with pytest.raises(RuntimeError) as excinfo:
        latest_ccle_names(["d", "dd", "ddd"], True, True)
    assert 'Mappings are not unique' in str(excinfo.value)
