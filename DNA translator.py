#!/usr/bin/env python3
"""
DNA to Protein Translator & Reverse Translator
Author: Ansuman Gadanayak
Description: A robust Python script to handle transcription, forward translation 
             (accounting for genetic code degeneracy), and reverse translation.
"""

# Standard Genetic Code Dictionary (Codon to Amino Acid)
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# Reverse mapping for degeneracy (Amino Acid -> List of Codons)
REVERSE_GENETIC_CODE = {}
for codon, aa in GENETIC_CODE.items():
    if aa not in REVERSE_GENETIC_CODE:
        REVERSE_GENETIC_CODE[aa] = []
    REVERSE_GENETIC_CODE[aa].append(codon)

def clean_dna(dna_sequence: str) -> str:
    """Clean whitespace and validate DNA nucleotides."""
    cleaned = "".join(dna_sequence.split()).upper()
    valid = set("ATCG")
    invalid = set(cleaned) - valid
    if invalid:
        raise ValueError(f"Invalid nucleotide(s): {', '.join(invalid)}. Use only A, T, C, G.")
    return cleaned

def clean_protein(protein_sequence: str) -> str:
    """Clean whitespace and validate protein amino acid letters."""
    cleaned = "".join(protein_sequence.split()).upper()
    valid_aas = set(GENETIC_CODE.values())
    invalid = set(cleaned) - valid_aas
    if invalid:
        raise ValueError(f"Invalid amino acid letter(s): {', '.join(invalid)}.")
    return cleaned

def translate_dna(dna_sequence: str, frame: int = 0) -> str:
    """Translate DNA to Protein across a specific reading frame (0, 1, or 2)."""
    cleaned = clean_dna(dna_sequence)
    if frame not in [0, 1, 2]:
        raise ValueError("Reading frame must be 0, 1, or 2.")
    
    protein = []
    shifted_seq = cleaned[frame:]
    
    for i in range(0, len(shifted_seq) - 2, 3):
        codon = shifted_seq[i:i+3]
        amino_acid = GENETIC_CODE.get(codon, 'X')
        protein.append(amino_acid)
        if amino_acid == '*':
            break
            
    return "".join(protein)

def reverse_translate_protein(protein_sequence: str) -> str:
    """
    Convert a Protein sequence back to a DNA sequence, 
    accounting for genetic code degeneracy by picking the first synonymous codon.
    """
    cleaned = clean_protein(protein_sequence)
    dna_pieces = []
    
    for aa in cleaned:
        if aa in REVERSE_GENETIC_CODE:
            dna_pieces.append(REVERSE_GENETIC_CODE[aa][0])
        else:
            dna_pieces.append("NNN")
            
    return "".join(dna_pieces)

if __name__ == "__main__":
    print("==================================================")
    print("   Bioinformatics DNA/Protein Translator Suite    ")
    print("==================================================\n")
    print("Instructions:")
    print("  - Paste DNA directly to translate (e.g., ATGGCCATT)")
    print("  - Type REV followed by protein sequence to reverse-translate (e.g., REV MAIV)")
    print("  - Type 'exit' or 'quit' to close.\n")
    
    while True:
        try:
            user_input = input("Enter sequence or command: ").strip()
            
            if not user_input:
                continue
            if user_input.lower() in ['exit', 'quit']:
                print("Exiting translator. Good luck with your GitHub repo!")
                break
                
            # Check if command starts with REV (case-insensitive)
            if user_input.upper().startswith("REV"):
                # Safely slice off 'REV' regardless of spacing
                parts = user_input.split(maxsplit=1)
                if len(parts) < 2:
                    print("\n[!] Error: Please provide a protein sequence after REV (e.g., REV MAIV)\n" + "-"*50 + "\n")
                    continue
                
                prot_seq = parts[1]
                dna_result = reverse_translate_protein(prot_seq)
                print(f"\n[+] Input Protein:  {clean_protein(prot_seq)}")
                print(f"[+] Degenerate DNA: {dna_result}\n" + "-"*50 + "\n")
                
            else:
                # Handle Forward Translation
                cleaned_dna = clean_dna(user_input)
                mrna = cleaned_dna.replace('T', 'U')
                protein = translate_dna(user_input, frame=0)
                
                print(f"\n[+] Length:     {len(cleaned_dna)} bp")
                print(f"[+] mRNA:       {mrna}")
                print(f"[+] Protein:    {protein} (Frame 0)")
                print("-" * 50 + "\n")
                
        except ValueError as e:
            print(f"\n[!] Error: {e}\n" + "-" * 50 + "\n")
        except KeyboardInterrupt:
            print("\nExiting translator. Goodbye!")
            break
