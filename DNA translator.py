
"""
DNA to Protein Translator and Sequence Analyzer
 By: Ansuman Gadanayak
Description: A Python tool featuring transcription, translation, reverse translation,
             GC content, melting temperature (Tm), molecular weight, ORFs, restriction mapping,
             primer characteristic checks, and codon optimization.
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

# High-frequency E. coli Codon Usage Table (for Codon Optimization)
E_COLI_OPTIMAL_CODONS = {
    'A': 'GCG', 'C': 'TGC', 'D': 'GAT', 'E': 'GAA', 'F': 'TTT',
    'G': 'GGC', 'H': 'CAT', 'I': 'ATT', 'K': 'AAA', 'L': 'CTG',
    'M': 'ATG', 'N': 'AAC', 'P': 'CCG', 'Q': 'CAG', 'R': 'CGT',
    'S': 'AGC', 'T': 'ACC', 'V': 'GTG', 'W': 'TGG', 'Y': 'TAT', '*': 'TAA'
}

# Common Restriction Enzymes and recognition sites
RESTRICTION_ENZYMES = {
    'EcoRI': 'GAATTC', 'BamHI': 'GGATCC', 'HindIII': 'AAGCTT',
    'NotI': 'GCGGCCGC', 'XbaI': 'TCTAGA', 'PstI': 'CTGCAG'
}


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


def calculate_gc_content(dna_sequence: str) -> float:
    """Calculate the GC content percentage."""
    cleaned = clean_dna(dna_sequence)
    if not cleaned:
        return 0.0
    return round(((cleaned.count('G') + cleaned.count('C')) / len(cleaned)) * 100, 2)


def calculate_melting_point(dna_sequence: str) -> float:
    """Calculate DNA melting temperature (Tm) via Wallace rule or empirical formula."""
    cleaned = clean_dna(dna_sequence)
    length = len(cleaned)
    if length == 0:
        return 0.0
    
    a, t, g, c = cleaned.count('A'), cleaned.count('T'), cleaned.count('G'), cleaned.count('C')
    if length <= 14:
        tm = (2 * (a + t)) + (4 * (g + c))
    else:
        tm = 64.9 + (41 * (g + c - 16.4) / length)
    return round(tm, 2)


def reverse_complement(dna_sequence: str) -> str:
    """Generate the reverse complement strand."""
    cleaned = clean_dna(dna_sequence)
    return cleaned.translate(str.maketrans("ATCG", "TAGC"))[::-1]


def calculate_molecular_weight(dna_sequence: str) -> float:
    """Calculate single-stranded DNA molecular weight in kDa."""
    cleaned = clean_dna(dna_sequence)
    if not cleaned:
        return 0.0
    weights = {'A': 313.21, 'T': 304.2, 'C': 289.18, 'G': 329.21}
    mw = sum(weights.get(nuc, 0) for nuc in cleaned) + 79.0
    return round(mw / 1000.0, 3)


def translate_dna(dna_sequence: str, frame: int = 0) -> str:
    """Translate DNA to protein for a selected reading frame."""
    cleaned = clean_dna(dna_sequence)
    if frame not in [0, 1, 2]:
        raise ValueError("Reading frame must be 0, 1, or 2.")
    protein = []
    for i in range(frame, len(cleaned) - 2, 3):
        codon = cleaned[i:i+3]
        aa = GENETIC_CODE.get(codon, 'X')
        protein.append(aa)
        if aa == '*':
            break
    return "".join(protein)


def reverse_translate_protein(protein_sequence: str) -> str:
    """Convert protein back to degenerate DNA using first choice codons."""
    cleaned = clean_protein(protein_sequence)
    return "".join(REVERSE_GENETIC_CODE[aa][0] if aa in REVERSE_GENETIC_CODE else "NNN" for aa in cleaned)


def optimize_for_ecoli(protein_sequence: str) -> str:
    """Optimize a protein sequence for high-level expression in E. coli."""
    cleaned = clean_protein(protein_sequence)
    return "".join(E_COLI_OPTIMAL_CODONS.get(aa, "NNN") for aa in cleaned)


def check_primer_properties(dna_sequence: str) -> dict:
    """Evaluate wet-lab PCR primer quality metrics (length, GC%, Tm, 3'-end clamp)."""
    cleaned = clean_dna(dna_sequence)
    length = len(cleaned)
    gc = calculate_gc_content(cleaned)
    tm = calculate_melting_point(cleaned)
    
    # Check 3' stability (primer should ideally end with a G or C to anchor well)
    ends_well = cleaned[-1] in {'G', 'C'} if length > 0 else False
    
    warnings = []
    if not (18 <= length <= 25):
        warnings.append("Length outside ideal PCR range (18-25 bp).")
    if not (40 <= gc <= 60):
        warnings.append("GC content outside optimal range (40-60%).")
    if not (52 <= tm <= 65):
        warnings.append("Tm outside optimal annealing window (52-65 °C).")
    if not ends_well:
        warnings.append("3'-end does not terminate in G or C (weak binding clamp).")

    return {
        'length': length,
        'gc': gc,
        'tm': tm,
        'ideal_clamp': ends_well,
        'warnings': warnings if warnings else ["None. Primer looks solid!"]
    }


def find_orfs(dna_sequence: str) -> list:
    """Scan all frames for open reading frames."""
    cleaned = clean_dna(dna_sequence)
    orfs = []
    for frame in [0, 1, 2]:
        shifted = cleaned[frame:]
        i = 0
        while i < len(shifted) - 2:
            if shifted[i:i+3] == 'ATG':
                start = frame + i
                j = i + 3
                codons = ['ATG']
                found_stop = False
                while j < len(shifted) - 2:
                    nxt = shifted[j:j+3]
                    codons.append(nxt)
                    if nxt in {'TAA', 'TAG', 'TGA'}:
                        found_stop = True
                        break
                    j += 3
                if found_stop:
                    seq = "".join(codons)
                    orfs.append({
                        'frame': f"+{frame+1}",
                        'start': start + 1,
                        'end': start + len(seq),
                        'length_bp': len(seq),
                        'protein': translate_dna(seq, 0)
                    })
                    i = j
            i += 3
    return orfs


def find_restriction_sites(dna_sequence: str) -> dict:
    """Locate restriction enzyme cutting positions."""
    cleaned = clean_dna(dna_sequence)
    results = {}
    for enzyme, site in RESTRICTION_ENZYMES.items():
        positions = []
        start = 0
        while True:
            pos = cleaned.find(site, start)
            if pos == -1:
                break
            positions.append(pos + 1)
            start = pos + 1
        if positions:
            results[enzyme] = positions
    return results


if __name__ == "__main__":
    print("==================================================")
    print("   Lab-Grade Bioinformatics DNA/Protein Suite     ")
    print("==================================================\n")
    print("Commands & Usage:")
    print("  - Paste DNA sequence directly for full sequence/primer analysis.")
    print("  - Type REV <protein> to reverse-translate.")
    print("  - Type OPT <protein> to codon-optimize for E. coli expression.")
    print("  - Type 'exit' or 'quit' to close.\n")

    while True:
        try:
            user_input = input("Enter sequence or command: ").strip()
            if not user_input:
                continue

            if user_input.lower() in ['exit', 'quit']:
                print("Exiting tool. Good luck with your lab experiments!")
                break

            cmd_upper = user_input.upper()

            if cmd_upper.startswith("REV"):
                parts = user_input.split(maxsplit=1)
                if len(parts) < 2:
                    print("\n[!] Error: Provide protein sequence (e.g., REV MAIV)\n" + "-"*50 + "\n")
                    continue
                dna_res = reverse_translate_protein(parts[1])
                print(f"\n[+] Input Protein:  {clean_protein(parts[1])}")
                print(f"[+] Degenerate DNA: {dna_res}\n" + "-"*50 + "\n")

            elif cmd_upper.startswith("OPT"):
                parts = user_input.split(maxsplit=1)
                if len(parts) < 2:
                    print("\n[!] Error: Provide protein sequence (e.g., OPT MAIV)\n" + "-"*50 + "\n")
                    continue
                opt_dna = optimize_for_ecoli(parts[1])
                print(f"\n[+] Input Protein:     {clean_protein(parts[1])}")
                print(f"[+] E. coli Opt DNA:   {opt_dna}\n" + "-"*50 + "\n")

            else:
                cleaned_dna = clean_dna(user_input)
                mrna = cleaned_dna.replace('T', 'U')
                protein = translate_dna(user_input, frame=0)
                gc = calculate_gc_content(cleaned_dna)
                tm = calculate_melting_point(cleaned_dna)
                rev_comp = reverse_complement(cleaned_dna)
                mw = calculate_molecular_weight(cleaned_dna)
                primer_check = check_primer_properties(cleaned_dna)
                orfs = find_orfs(cleaned_dna)
                restriction_sites = find_restriction_sites(cleaned_dna)

                print(f"\n[+] Length:          {len(cleaned_dna)} bp")
                print(f"[+] mRNA:            {mrna}")
                print(f"[+] Rev Complement:  {rev_comp}")
                print(f"[+] Protein (F0):    {protein}")
                print(f"[+] GC Content:      {gc}%")
                print(f"[+] Melting Pt (Tm): {tm} °C")
                print(f"[+] Molecular Wt:    {mw} kDa")
                
                print("\n[+] PCR Primer Quality Metrics:")
                print(f"    - Ideal 3' Clamp: {'Yes (G/C)' if primer_check['ideal_clamp'] else 'No (A/T)'}")
                print("    - Evaluation Notes:")
                for warning in primer_check['warnings']:
                    print(f"      * {warning}")

                print("\n[+] Open Reading Frames (ORFs):")
                if orfs:
                    for idx, orf in enumerate(orfs, 1):
                        print(f"    {idx}. Frame {orf['frame']} | Pos: {orf['start']}-{orf['end']} ({orf['length_bp']} bp)")
                else:
                    print("    None detected.")

                print("\n[+] Restriction Enzyme Sites:")
                if restriction_sites:
                    for enzyme, pos_list in restriction_sites.items():
                        print(f"    - {enzyme}: {pos_list}")
                else:
                    print("    None of the configured restriction sites found.")

                print("-" * 50 + "\n")

        except ValueError as e:
            print(f"\n[!] Error: {e}\n" + "-" * 50 + "\n")
        except KeyboardInterrupt:
            print("\nExiting tool. Goodbye!")
            break
    
    
    
    
       
               
              
