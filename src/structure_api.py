import requests

def fetch_pdb_file(pdb_id: str) -> str | None:
    """Fetches PDB file from RCSB database."""
    url = f"https://files.rcsb.org/download/{pdb_id.lower()}.pdb"
    res = requests.get(url)
    return res.text if res.status_code == 200 else None

def fetch_string_interactions(protein_name: str) -> list:
    """Queries STRING DB for functional interaction partners."""
    url = "https://string-db.org/api/json/network"
    params = {"identifiers": protein_name, "species": 9606, "limit": 5}
    res = requests.get(url, params=params)
    return res.json() if res.status_code == 200 else []
