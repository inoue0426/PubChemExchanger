import argparse
import pandas as pd
import requests
from functools import lru_cache

@lru_cache(maxsize=None)
def fetch_data(url):
    response = requests.get(url)
    return response.json() if response.status_code == 200 else None

def get_property(data, key, default="Not found"):
    return data.get(key, default)

def get_identifiers(data, key):
    return (
        data["IdentifierList"].get(key, []) if data and "IdentifierList" in data else []
    )

def get_synonyms(data):
    return (
        data["InformationList"]["Information"][0].get("Synonym", [])
        if data
        else []
    )

def extract_nsc_numbers(synonyms):
    return list(set("".join(filter(str.isdigit, syn)) for syn in synonyms if syn.startswith("NSC")))

def extract_dbid(synonyms):
    return next((syn for syn in synonyms if syn.startswith("DB")), "Not found")

def get_drug_info(drug_names):
    drugs_data = []
    for drug_name in drug_names:
        data = {"SMILES": "Not found", "CID": [], "SID": [], "NSC": [], "DBID": "Not found", "Synonyms": []}
        base_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/"

        properties = fetch_data(f"{base_url}property/CanonicalSMILES/JSON")
        data["SMILES"] = get_property(properties["PropertyTable"]["Properties"][0], "CanonicalSMILES") if properties else "Not found"

        data["CID"] = get_identifiers(fetch_data(f"{base_url}cids/JSON"), "CID")
        data["SID"] = get_identifiers(fetch_data(f"{base_url}sids/JSON"), "SID")

        synonyms = fetch_data(f"{base_url}synonyms/JSON")
        data["Synonyms"] = get_synonyms(synonyms)
        data["NSC"] = extract_nsc_numbers(data["Synonyms"])
        data["DBID"] = extract_dbid(data["Synonyms"])

        drugs_data.append(data)

    return pd.DataFrame(drugs_data, index=drug_names)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve information about drugs from PubChem")
    parser.add_argument("-d", "--drug_names", nargs='+', type=str, help="The names of the drugs to query")
    args = parser.parse_args()
    result = get_drug_info(args.drug_names)
    print(result)
