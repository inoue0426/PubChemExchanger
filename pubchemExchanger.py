import argparse

import pandas as pd
import requests


def fetch_data(url):
    response = requests.get(url)
    return response.json() if response.status_code == 200 else None


def get_property(data, key, default="Not found"):
    return data.get(key, default)


def get_identifiers(data, key):
    return (
        data["IdentifierList"].get(key, []) if data and "IdentifierList" in data else []
    )


def get_synonyms(data, prefix):
    return (
        [
            syn
            for syn in data["InformationList"]["Information"][0]["Synonym"]
            if syn.startswith(prefix)
        ]
        if data
        else []
    )


def get_drug_info(drug_name):
    data = {"SMILES": "Not found", "CID": [], "SID": [], "NSC": [], "DBID": "Not found"}

    base_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/"

    properties = fetch_data(f"{base_url}property/CanonicalSMILES/JSON")
    data["SMILES"] = (
        get_property(properties["PropertyTable"]["Properties"][0], "CanonicalSMILES")
        if properties
        else "Not found"
    )

    data["CID"] = get_identifiers(fetch_data(f"{base_url}cids/JSON"), "CID")
    data["SID"] = get_identifiers(fetch_data(f"{base_url}sids/JSON"), "SID")

    synonyms = fetch_data(f"{base_url}synonyms/JSON")
    nsc_synonyms = get_synonyms(synonyms, "NSC")
    data["NSC"] = list(set("".join(filter(str.isdigit, syn)) for syn in nsc_synonyms))
    data["DBID"] = next((syn for syn in get_synonyms(synonyms, "DB")), "Not found")

    return pd.DataFrame([data], index=[drug_name])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Retrieve information about a drug from PubChem"
    )
    parser.add_argument(
        "-d", "--drug_name", type=str, help="The name of the drug to query"
    )
    args = parser.parse_args()
    result = get_drug_info(args.drug_name)
    print(result)
