
# DNA-TO-PROTEIN 🧬

A simple, robust Python command-line tool for bioinformatics that handles DNA transcription, forward translation into proteins, and reverse translation back to DNA.

## Features
 Forward Translation: Converts DNA sequences into protein sequences across reading frames.
 Reverse Translation: Converts protein sequences back into degenerate DNA sequences.
 Degeneracy Handling: Accounts for the fact that multiple codons can code for the same amino acid.
 Interactive Shell: Keeps the terminal open so you can test multiple sequences back-to-back without restarting.

---

## Understanding the Output & Symbols

When you run translations, you might notice a few special characters:
 * (Stop Codon):* Represents a stop signal (like UAA, UAG, or UGA) which tells the cell to stop translating the protein.
X (Unknown/Invalid):** Used as a fallback if an unrecognized character or codon slips through.

---

## Why Genetic Code Degeneracy Matters (And What i Did)

In biology, the genetic code is degenerate, meaning most amino acids can be encoded by more than one codon (for example, Leucine can be written as six different codons like CTT, CTC, CTA, CTG, TTA, or TTG). 

 The Problem: When you go backwards from a protein to DNA (reverse translation), a standard program doesn't know *which* original codon to pick.
 Our Solution: We built a custom REVERSE GENETIC CODE dictionary that maps every single amino acid to its valid synonymous codons, allowing the script to cleanly map proteins back into functional DNA sequences.

---
