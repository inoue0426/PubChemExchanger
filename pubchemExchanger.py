import argparse
import time
from functools import lru_cache

import pandas as pd
import requests
from joblib import Parallel, delayed
from tqdm import tqdm


@lru_cache(maxsize=None)
def fetch_data(url):
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Request failed: {e}, retrying.")
        time.sleep(2)
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Retry failed: {e}")
            return None


def _get_property(data, key, default="Not found"):
    return data.get(key, default)


def _get_identifiers(data, key):
    return (
        data["IdentifierList"].get(key, []) if data and "IdentifierList" in data else []
    )


def _get_synonyms(data):
    return data["InformationList"]["Information"][0].get("Synonym", []) if data else []


def _extract_nsc_numbers(synonyms):
    return list(
        set(
            "".join(filter(str.isdigit, syn))
            for syn in synonyms
            if syn.startswith("NSC")
        )
    )


def _extract_dbid(synonyms):
    for syn in synonyms:
        if syn.startswith("DB"):
            return syn.replace("-", "")
    return "Not found"


def process_drug(identifier, identifier_type="name"):
    data = {
        "SMILES": pd.NA,
        "CID": [],
        "SID": [],
        "NSC": [],
        "DBID": pd.NA,
        "Synonyms": [],
        "ChEMBL_ID": pd.NA,
        "PubChem_ID": pd.NA,
    }
    base_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{identifier_type}/{identifier}/"

    properties = fetch_data(
        f"{base_url}property/CanonicalSMiles,ChEMBL_ID,PubChem_ID/JSON"
    )
    if properties:
        properties_data = properties["PropertyTable"]["Properties"][0]
        data["SMILES"] = _get_property(properties_data, "CanonicalSMILES")
        data["ChEMBL_ID"] = _get_property(properties_data, "ChEMBL_ID")
        data["PubChem_ID"] = _get_property(properties_data, "PubChem_ID")

    data["CID"] = _get_identifiers(fetch_data(f"{base_url}cids/JSON"), "CID")
    data["SID"] = _get_identifiers(fetch_data(f"{base_url}sids/JSON"), "SID")

    synonyms = fetch_data(f"{base_url}synonyms/JSON")
    data["Synonyms"] = _get_synonyms(synonyms)
    data["NSC"] = _extract_nsc_numbers(data["Synonyms"])
    data["DBID"] = _extract_dbid(data["Synonyms"])

    return data


def get_drug_info(identifiers, identifier_type="name"):
    drugs_data = Parallel(n_jobs=-1)(
        delayed(process_drug)(identifier, identifier_type)
        for identifier in tqdm(identifiers, desc="Processing drugs")
    )
    return pd.DataFrame(drugs_data, index=identifiers)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Retrieve information about drugs from PubChem"
    )
    parser.add_argument(
        "-d",
        "--drug_names",
        nargs="+",
        type=str,
        help="The names of the drugs to query",
    )
    parser.add_argument(
        "-t",
        "--type",
        default="name",
        choices=["name", "cid", "sid", "chembl_id"],
        help="Type of identifier",
    )
    args = parser.parse_args()
    result = get_drug_info(args.drug_names, args.type)
    print(result)
