# DNA-TO-PROTEIN 
<img width="1920" height="1019" alt="Screenshot (13)" src="https://github.com/user-attachments/assets/f325634e-7d77-45ea-913d-c80840fd048b" />

 Python command line tool for bioinformatics that handles DNA transcription, forward translation into proteins, and reverse translation back to DNA.

## Features

* Forward Translation: Converts DNA sequences into protein sequences across reading frames.
  
* Reverse Translation: Converts protein sequences back into degenerate DNA sequences.
  
* Degeneracy Handling: Accounts for the fact that multiple codons can code for the same amino acid.
  
* Interactive Shell: Keeps the terminal open so you can test multiple sequences back to back without restarting.
  
* GC Content Calculation: Calculates the GC percentage to measure sequence stability.
  
* Melting Point (Tm) Estimation: Computes the melting temperature using the Wallace rule for short sequences and the empirical formula for longer oligos.
  
* Reverse Complement Generator: Generates the complementary antiparallel DNA strand automatically.
  
* Molecular Weight: Computes the approximate molecular weight of single stranded DNA in kDa.
  
* PCR Primer QC Checks: Checks length, GC percentage, Tm range, and 3'-end anchoring to evaluate primer design quality.
 
* Open Reading Frame (ORF) Finder: Scans reading frames for start codons (ATG) and in-frame stop codons.
 
* Restriction Site Mapping: Finds cut positions for common restriction enzymes (EcoRI, BamHI, HindIII, NotI, XbaI, PstI).
  
* Codon Optimization for E. coli: Converts amino acid sequences into preferred bacterial codons for recombinant protein expression.

## Understanding the Output and Symbols

When you run translations, you might notice a few special characters:
 * (Stop Codon):* Represents a stop signal (like TAA, TAG, or TGA) which tells the cell to stop translating the protein.

* X (Unknown/Invalid): Used as a fallback value if an unknown character or wrong codon gets through.
  
* NNN: Returned during reverse translation if an unknown amino acid is provided.

## Why Genetic Code Degeneracy Matters (And What I Did)

 the genetic code is degenerate, meaning most amino acids can be encoded by more than one codon (for example, Leucine can be written as six different codons like CTT, CTC, CTA, CTG, TTA, or TTG).

*  Problem: When you go backwards from a protein to DNA (reverse translation), a standard program doesn't know *which* original codon to pick.
 
*  Solution: made a custom reverse genetic code dictionary that maps every single amino acid to its valid synonymous codons, allowing the script to cleanly map proteins back into functional DNA sequences.
