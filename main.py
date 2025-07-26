#!/usr/bin/env python3
"""
Protein Enrichment Analysis Script

This script analyzes the enrichment of amino acids at the 2nd and 3rd positions
in proteins identified in a screen compared to the human proteome.

The script works with:
- Excel file containing gene names identified in the screen (files/UBR3PEP.xlsx)
- FASTA file containing human proteome sequences (files/PROTEOME.fasta)

Requirements:
- biopython: pip install biopython
- scipy: pip install scipy
- pandas: pip install pandas
- statsmodels: pip install statsmodels
- openpyxl: pip install openpyxl

Usage:
1. Create a 'files/' folder and place your input files there:
   - files/UBR3PEP.xlsx (your screen results)
   - files/PROTEOME.fasta (human proteome from UniProt)
2. Run: python peptide_enrichment_analysis.py
3. Results will be saved in the 'outputs/' folder
"""

from Bio import SeqIO
from collections import Counter
from scipy.stats import fisher_exact
import pandas as pd
from statsmodels.stats.multitest import multipletests
import os
import sys
import re


def read_fasta_sequences(fasta_file):
    """
    Read protein sequences from a FASTA file and create a mapping of gene names to sequences.

    Args:
        fasta_file (str): Path to the FASTA file

    Returns:
        tuple: (all_sequences, gene_to_sequence_map)
    """
    all_sequences = []
    gene_to_sequence = {}

    try:
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequence = str(record.seq)
                all_sequences.append(sequence)

                # Extract gene name from the header
                header = record.description
                gene_name = extract_gene_name(header)
                if gene_name:
                    # Store the sequence (if multiple isoforms, keep the first one)
                    if gene_name not in gene_to_sequence:
                        gene_to_sequence[gene_name] = sequence

        print(f"Successfully read {len(all_sequences)} protein sequences from {fasta_file}")
        print(f"Mapped {len(gene_to_sequence)} unique gene names to sequences")
        return all_sequences, gene_to_sequence

    except FileNotFoundError:
        print(f"Error: Could not find file {fasta_file}")
        return [], {}
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return [], {}


def extract_gene_name(header):
    """
    Extract gene name from UniProt FASTA header.

    Args:
        header (str): FASTA header line

    Returns:
        str: Gene name or None
    """
    # Try to extract gene name from GN= field
    gn_match = re.search(r'GN=([^\s]+)', header)
    if gn_match:
        return gn_match.group(1)

    # If no GN field, try to extract from the protein name part
    # Format: >sp|ID|NAME_HUMAN or >tr|ID|NAME_HUMAN
    parts = header.split('|')
    if len(parts) >= 3:
        name_part = parts[2].split('_')[0]  # Take part before _HUMAN
        return name_part

    return None


def read_screen_genes(excel_file):
    """
    Read gene names from the Excel file.

    Args:
        excel_file (str): Path to the Excel file

    Returns:
        list: List of gene names identified in the screen
    """
    try:
        # Read the Excel file
        df = pd.read_excel(excel_file)

        print(f"Excel file columns: {list(df.columns)}")
        print(f"First few rows:")
        print(df.head())

        # The gene names appear to be in the third column (index 2)
        # Based on the data structure we saw: [condition, transcript_id, gene_name, pdf_file]
        gene_column = df.iloc[:, 2]  # Third column (0-indexed)
        gene_names = gene_column.dropna().unique().tolist()

        print(f"\nExtracted {len(gene_names)} unique gene names from screen")
        print(f"Sample gene names: {gene_names[:10]}")

        return gene_names

    except FileNotFoundError:
        print(f"Error: Could not find file {excel_file}")
        return []
    except Exception as e:
        print(f"Error reading Excel file: {e}")
        return []


def get_screen_sequences(screen_genes, gene_to_sequence):
    """
    Get protein sequences for genes identified in the screen.

    Args:
        screen_genes (list): List of gene names from screen
        gene_to_sequence (dict): Mapping of gene names to sequences

    Returns:
        tuple: (found_sequences, not_found_genes)
    """
    found_sequences = []
    not_found_genes = []

    for gene in screen_genes:
        if gene in gene_to_sequence:
            found_sequences.append(gene_to_sequence[gene])
        else:
            not_found_genes.append(gene)

    print(f"\nSequence matching results:")
    print(f"Found sequences for {len(found_sequences)} genes")
    print(f"Could not find sequences for {len(not_found_genes)} genes")

    if not_found_genes:
        print(f"Genes not found: {not_found_genes[:10]}{'...' if len(not_found_genes) > 10 else ''}")

    return found_sequences, not_found_genes


def analyze_residue_frequencies(sequences, sequence_type="sequences"):
    """
    Analyze amino acid frequencies at 2nd and 3rd positions.

    Args:
        sequences (list): List of protein sequences
        sequence_type (str): Description for output messages

    Returns:
        tuple: (second_residue_counts, third_residue_counts, total_with_2nd, total_with_3rd)
    """
    second_residue_counts = Counter()
    third_residue_counts = Counter()
    total_with_second = 0
    total_with_third = 0

    print(f"\nAnalyzing {len(sequences)} {sequence_type}...")

    for seq in sequences:
        # Check for second residue (position 1 in 0-indexed)
        if len(seq) >= 2:
            second_residue = seq[1]
            if second_residue.isalpha():  # Make sure it's a valid amino acid
                second_residue_counts[second_residue] += 1
                total_with_second += 1

        # Check for third residue (position 2 in 0-indexed)
        if len(seq) >= 3:
            third_residue = seq[2]
            if third_residue.isalpha():  # Make sure it's a valid amino acid
                third_residue_counts[third_residue] += 1
                total_with_third += 1

    print(f"Total {sequence_type} with 2nd residue: {total_with_second}")
    print(f"Total {sequence_type} with 3rd residue: {total_with_third}")

    return second_residue_counts, third_residue_counts, total_with_second, total_with_third


def print_frequency_summary(residue_counts, total_count, position_name, dataset_name):
    """Print a summary of amino acid frequencies."""
    print(f"\n{dataset_name} - {position_name} Residue Frequencies:")
    print("-" * 60)
    for aa in sorted(residue_counts.keys()):
        count = residue_counts[aa]
        frequency = count / total_count if total_count > 0 else 0
        print(f"{aa}: {count:>6} ({frequency:.4f})")


def perform_enrichment_analysis(screen_counts, proteome_counts, total_screen, total_proteome, position_name):
    """
    Perform enrichment analysis for a specific position.

    Args:
        screen_counts (Counter): Amino acid counts from screen
        proteome_counts (Counter): Amino acid counts from proteome
        total_screen (int): Total sequences with this position in screen
        total_proteome (int): Total sequences with this position in proteome
        position_name (str): Name of the position being analyzed

    Returns:
        pandas.DataFrame: Results of enrichment analysis
    """
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"  # 20 standard amino acids
    results = []
    p_values = []

    print(f"\n--- Enrichment Analysis for {position_name} Residue ---")

    for aa in amino_acids:
        count_screen = screen_counts.get(aa, 0)
        count_proteome = proteome_counts.get(aa, 0)

        if total_screen == 0 or total_proteome == 0:
            fold_enrichment = float('nan')
            p_value = float('nan')
        else:
            freq_screen = count_screen / total_screen
            freq_proteome = count_proteome / total_proteome

            # Calculate fold enrichment
            if freq_proteome == 0:
                fold_enrichment = float('inf') if freq_screen > 0 else 1.0
            else:
                fold_enrichment = freq_screen / freq_proteome

            # Fisher's Exact Test
            table = [
                [count_screen, total_screen - count_screen],
                [count_proteome, total_proteome - count_proteome]
            ]

            try:
                odds_ratio, p_value = fisher_exact(table, alternative='two-sided')
            except Exception as e:
                print(f"Error in Fisher's exact test for {aa}: {e}")
                p_value = float('nan')

        results.append({
            'Amino_Acid': aa,
            'Count_Screen': count_screen,
            'Freq_Screen': freq_screen if total_screen > 0 else 0,
            'Count_Proteome': count_proteome,
            'Freq_Proteome': freq_proteome if total_proteome > 0 else 0,
            'Fold_Enrichment': fold_enrichment,
            'P_Value': p_value
        })
        p_values.append(p_value)

    # Convert to DataFrame
    df = pd.DataFrame(results)

    # Multiple testing correction (Benjamini-Hochberg)
    valid_p_values = [p for p in p_values if pd.notna(p)]
    if valid_p_values:
        try:
            reject, pvals_corrected, _, _ = multipletests(valid_p_values, method='fdr_bh')
            corrected_p_map = {orig_p: corr_p for orig_p, corr_p in zip(valid_p_values, pvals_corrected)}
            df['FDR_P_Value'] = df['P_Value'].apply(lambda p: corrected_p_map.get(p, float('nan')))
            df['Significant'] = df['FDR_P_Value'] < 0.05
        except Exception as e:
            print(f"Error in multiple testing correction: {e}")
            df['FDR_P_Value'] = float('nan')
            df['Significant'] = False
    else:
        df['FDR_P_Value'] = float('nan')
        df['Significant'] = False

    return df


def save_results_to_csv(df_second, df_third, screen_genes, not_found_genes,
                        output_file="outputs/UBR3_enrichment_results.csv"):
    """Save results to CSV file with additional metadata."""
    try:
        # Create outputs directory if it doesn't exist
        os.makedirs("outputs", exist_ok=True)

        # Add position column to distinguish results
        df_second_copy = df_second.copy()
        df_second_copy['Position'] = 'Second'
        df_third_copy = df_third.copy()
        df_third_copy['Position'] = 'Third'

        # Combine results
        combined_df = pd.concat([df_second_copy, df_third_copy], ignore_index=True)
        combined_df.to_csv(output_file, index=False)

        # Save metadata
        metadata_file = output_file.replace('.csv', '_metadata.txt')
        with open(metadata_file, 'w') as f:
            f.write("UBR3 Peptide Screen Enrichment Analysis\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Total genes in screen: {len(screen_genes)}\n")
            f.write(f"Genes with sequences found: {len(screen_genes) - len(not_found_genes)}\n")
            f.write(f"Genes not found in proteome: {len(not_found_genes)}\n\n")

            if not_found_genes:
                f.write("Genes not found in proteome:\n")
                for gene in not_found_genes:
                    f.write(f"- {gene}\n")

        # Save simplified enrichment summary
        enrichment_summary_file = "outputs/UBR3_enrichment_summary.csv"
        summary_data = []

        for _, row in combined_df.iterrows():
            summary_data.append({
                'Position': row['Position'],
                'Amino_Acid': row['Amino_Acid'],
                'Screen_Frequency': row['Freq_Screen'],
                'Proteome_Frequency': row['Freq_Proteome'],
                'Enrichment_Ratio': row['Fold_Enrichment'],
                'Is_Enriched': 'Yes' if row['Fold_Enrichment'] > 1 else 'No',
                'Significant': 'Yes' if row['Significant'] else 'No',
                'FDR_P_Value': row['FDR_P_Value']
            })

        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(enrichment_summary_file, index=False)

        print(f"\nResults saved to {output_file}")
        print(f"Metadata saved to {metadata_file}")
        print(f"Enrichment summary saved to {enrichment_summary_file}")

    except Exception as e:
        print(f"Error saving results: {e}")


def main():
    """Main function to run the enrichment analysis."""
    print("=" * 80)
    print("UBR3 PEPTIDE SCREEN ENRICHMENT ANALYSIS")
    print("=" * 80)

    # File paths
    human_proteome_file = "files/PROTEOME.fasta"
    screen_excel_file = "files/UBR3PEP.xlsx"

    # Check if files exist
    if not os.path.exists(human_proteome_file):
        print(f"Error: Human proteome file '{human_proteome_file}' not found!")
        print("Please make sure your FASTA file is named 'PROTEOME.fasta' in the 'files/' folder")
        return

    if not os.path.exists(screen_excel_file):
        print(f"Error: Screen Excel file '{screen_excel_file}' not found!")
        print("Please make sure your Excel file is named 'UBR3PEP.xlsx' in the 'files/' folder")
        return

    # Step 1: Read the human proteome
    print("Step 1: Reading human proteome FASTA file...")
    all_proteome_sequences, gene_to_sequence = read_fasta_sequences(human_proteome_file)

    if not all_proteome_sequences:
        print("Error: Could not read proteome sequences.")
        return

    # Step 2: Read the screen genes from Excel
    print("\nStep 2: Reading screen genes from Excel file...")
    screen_genes = read_screen_genes(screen_excel_file)

    if not screen_genes:
        print("Error: Could not read screen genes.")
        return

    # Step 3: Get sequences for screen genes
    print("\nStep 3: Matching screen genes to protein sequences...")
    screen_sequences, not_found_genes = get_screen_sequences(screen_genes, gene_to_sequence)

    if not screen_sequences:
        print("Error: No matching sequences found for screen genes.")
        return

    # Step 4: Analyze proteome
    print("\nStep 4: Analyzing human proteome...")
    proteome_2nd, proteome_3rd, proteome_total_2nd, proteome_total_3rd = analyze_residue_frequencies(
        all_proteome_sequences, "proteome proteins"
    )

    # Step 5: Analyze screen proteins
    print("\nStep 5: Analyzing screen proteins...")
    screen_2nd, screen_3rd, screen_total_2nd, screen_total_3rd = analyze_residue_frequencies(
        screen_sequences, "screen proteins"
    )

    # Print frequency summaries
    print_frequency_summary(proteome_2nd, proteome_total_2nd, "Second", "Human Proteome")
    print_frequency_summary(screen_2nd, screen_total_2nd, "Second", "UBR3 Screen")
    print_frequency_summary(proteome_3rd, proteome_total_3rd, "Third", "Human Proteome")
    print_frequency_summary(screen_3rd, screen_total_3rd, "Third", "UBR3 Screen")

    # Step 6: Perform enrichment analysis
    print("\nStep 6: Performing enrichment analysis...")

    # Analyze second residue
    results_second = perform_enrichment_analysis(
        screen_2nd, proteome_2nd, screen_total_2nd, proteome_total_2nd, "Second"
    )

    # Analyze third residue
    results_third = perform_enrichment_analysis(
        screen_3rd, proteome_3rd, screen_total_3rd, proteome_total_3rd, "Third"
    )

    # Display results
    print("\n" + "=" * 80)
    print("ENRICHMENT RATIOS - SCREEN FREQUENCY ÷ PROTEOME FREQUENCY")
    print("=" * 80)
    print("This shows how many times more/less frequent each amino acid is")
    print("in your screen compared to the human proteome background")
    print("-" * 80)

    # Create a summary table for easy interpretation
    enrichment_summary = []

    for _, row in results_second.iterrows():
        enrichment_summary.append({
            'Position': 'Second',
            'Amino_Acid': row['Amino_Acid'],
            'Screen_Freq': f"{row['Freq_Screen']:.4f}",
            'Proteome_Freq': f"{row['Freq_Proteome']:.4f}",
            'Enrichment_Ratio': f"{row['Fold_Enrichment']:.2f}",
            'Interpretation': 'Enriched' if row['Fold_Enrichment'] > 1 else 'Depleted',
            'Significant': 'Yes' if row['Significant'] else 'No'
        })

    for _, row in results_third.iterrows():
        enrichment_summary.append({
            'Position': 'Third',
            'Amino_Acid': row['Amino_Acid'],
            'Screen_Freq': f"{row['Freq_Screen']:.4f}",
            'Proteome_Freq': f"{row['Freq_Proteome']:.4f}",
            'Enrichment_Ratio': f"{row['Fold_Enrichment']:.2f}",
            'Interpretation': 'Enriched' if row['Fold_Enrichment'] > 1 else 'Depleted',
            'Significant': 'Yes' if row['Significant'] else 'No'
        })

    summary_df = pd.DataFrame(enrichment_summary)

    # Show second position enrichments
    second_summary = summary_df[summary_df['Position'] == 'Second'].copy()
    second_summary = second_summary.sort_values('Enrichment_Ratio', ascending=False,
                                                key=lambda x: pd.to_numeric(x, errors='coerce'))
    print("\nSECOND POSITION ENRICHMENT RATIOS:")
    print(second_summary.to_string(index=False))

    # Show third position enrichments
    third_summary = summary_df[summary_df['Position'] == 'Third'].copy()
    third_summary = third_summary.sort_values('Enrichment_Ratio', ascending=False,
                                              key=lambda x: pd.to_numeric(x, errors='coerce'))
    print("\nTHIRD POSITION ENRICHMENT RATIOS:")
    print(third_summary.to_string(index=False))

    print("\n" + "=" * 80)
    print("DETAILED STATISTICAL RESULTS")
    print("=" * 80)
    print("ENRICHMENT RESULTS - SECOND RESIDUE")
    print("=" * 80)
    print("Sorted by FDR-corrected P-value (most significant first)")
    print("-" * 80)
    second_sorted = results_second.sort_values('FDR_P_Value').reset_index(drop=True)
    print(second_sorted.to_string(index=False, float_format='%.4f'))

    print("\n" + "=" * 80)
    print("ENRICHMENT RESULTS - THIRD RESIDUE")
    print("=" * 80)
    print("Sorted by FDR-corrected P-value (most significant first)")
    print("-" * 80)
    third_sorted = results_third.sort_values('FDR_P_Value').reset_index(drop=True)
    print(third_sorted.to_string(index=False, float_format='%.4f'))

    # Highlight significant results
    print("\n" + "=" * 80)
    print("SIGNIFICANT ENRICHMENTS (FDR p-value < 0.05)")
    print("=" * 80)

    sig_second = second_sorted[second_sorted['Significant'] == True]
    sig_third = third_sorted[third_sorted['Significant'] == True]

    if len(sig_second) > 0:
        print("\nSignificant Second Residue Enrichments:")
        print(sig_second[['Amino_Acid', 'Fold_Enrichment', 'FDR_P_Value']].to_string(index=False, float_format='%.4f'))
    else:
        print("\nNo significant enrichments found for second residue.")

    if len(sig_third) > 0:
        print("\nSignificant Third Residue Enrichments:")
        print(sig_third[['Amino_Acid', 'Fold_Enrichment', 'FDR_P_Value']].to_string(index=False, float_format='%.4f'))
    else:
        print("\nNo significant enrichments found for third residue.")

    # Save results
    print("\nStep 7: Saving results...")
    save_results_to_csv(results_second, results_third, screen_genes, not_found_genes)

    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE!")
    print("=" * 80)
    print(f"Analyzed {len(screen_sequences)} proteins from your UBR3 screen")
    print(f"against {len(all_proteome_sequences)} proteins in the human proteome")
    print("\nInterpretation guide:")
    print("- Fold_Enrichment > 1: Amino acid is more frequent in your screen than in proteome")
    print("- Fold_Enrichment < 1: Amino acid is less frequent in your screen than in proteome")
    print("- FDR_P_Value < 0.05: Statistically significant result")
    print("- Check the 'Significant' column for amino acids with reliable enrichment/depletion")


if __name__ == "__main__":
    main()