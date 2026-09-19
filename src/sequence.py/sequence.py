from Bio.SeqUtils.ProtParam import ProteinAnalysis

def parse_fasta(raw_text: str) -> tuple[str, str]:
    """Parses input text to separate FASTA header from sequence."""
    lines = raw_text.strip().splitlines()
    if not lines:
        return "Unknown Protein", ""
    
    if lines[0].startswith(">"):
        header = lines[0][1:].strip()
        sequence = "".join(lines[1:]).replace(" ", "").upper()
    else:
        header = "User Provided Sequence"
        sequence = "".join(lines).replace(" ", "").upper()
        
    return header, sequence

def validate_sequence(sequence: str) -> tuple[bool, str]:
    """Validates standard 20 amino acids."""
    if not sequence:
        return False, "Sequence is empty."
    
    valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")
    invalid_chars = set(sequence.upper()) - valid_amino_acids
    
    if invalid_chars:
        return False, f"Invalid amino acids: {', '.join(invalid_chars)}"
    
    return True, "Valid amino acid sequence."

def calculate_sequence_properties(sequence: str) -> dict:
    """Calculates MW, pI, GRAVY, composition, and charge counts."""
    analysis = ProteinAnalysis(sequence)
    
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    raw_counts = analysis.count_amino_acids()
    seq_len = len(sequence)
    composition = {aa: (raw_counts.get(aa, 0) / seq_len) * 100 for aa in amino_acids}
    
    positively_charged = raw_counts.get("K", 0) + raw_counts.get("R", 0) + raw_counts.get("H", 0)
    negatively_charged = raw_counts.get("D", 0) + raw_counts.get("E", 0)
    neutral = seq_len - (positively_charged + negatively_charged)
    
    return {
        "length": seq_len,
        "molecular_weight": round(analysis.molecular_weight(), 2),
        "isoelectric_point": round(analysis.isoelectric_point(), 2),
        "gravy": round(analysis.gravy(), 2),
        "composition": composition,
        "charge_distribution": {
            "Positive (+): K, R, H": positively_charged,
            "Negative (-): D, E": negatively_charged,
            "Neutral": neutral
        }
    }