# gRNA_design.py

from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import subprocess
import os
import re

# Define the path to your Bowtie index files
# This dictionary must be in this file.
GENOME_INDEXES = {
    "human": "/path/to/your/human_hg38_index",
    "mouse": "/path/to/your/mouse_mm10_index",
    # Add more organisms here
}
def run_bowtie_search(spacer_sequence, genome_index_path, max_mismatches=3):
    """
    Runs a Bowtie search to find off-target sites with mismatch details.
    
    Returns:
        list: A list of mismatch counts for each off-target site found.
    """
    clean_spacer = spacer_sequence.replace('N', 'A')
    
    command = [
        "bowtie",
        "-n", str(max_mismatches),
        "-l", "10",
        "--norc",
        "-a",
        "--sam", # Output in SAM format
        genome_index_path,
        "-",
        os.devnull
    ]

    mismatch_counts = []
    
    try:
        p = subprocess.run(
            command,
            input=f">spacer\n{clean_spacer}\n".encode('utf-8'),
            check=True,
            capture_output=True,
            text=True,
            timeout=60
        )
        
        for line in p.stdout.splitlines():
            if line.startswith('@'): continue
            
            sam_fields = line.split('\t')
            optional_fields = sam_fields[11:]
            
            mismatch_pattern = re.compile(r'NM:i:(\d+)')
            
            for field in optional_fields:
                match = mismatch_pattern.match(field)
                if match:
                    mismatch_counts.append(int(match.group(1)))
                    break
        
        return mismatch_counts

    except subprocess.CalledProcessError as e:
        print(f"Error running Bowtie: {e.stderr}")
        return []
    except FileNotFoundError:
        print("Error: Bowtie command not found. Is it in your PATH?")
        return []
    except subprocess.TimeoutExpired:
        print("Bowtie command timed out.")
        return []

def find_gRNAs(dna_sequence, pam_sequence="NGG", spacer_length=20):
    """
    Scans a DNA sequence for all potential gRNA target sites based on a PAM.
    
    
    """
    gRNAs = []
    
    # Biopython's Seq object handles DNA and its reverse complement easily
    dna_seq = Seq(dna_sequence.upper())
    
    # 1. Search for gRNAs on the forward strand
    for i in range(len(dna_seq) - spacer_length - len(pam_sequence) + 1):
        potential_pam = dna_seq[i + spacer_length : i + spacer_length + len(pam_sequence)]
        is_pam = (
            (pam_sequence[0] == 'N' or potential_pam[0] == pam_sequence[0]) and
            (pam_sequence[1] == 'N' or potential_pam[1] == pam_sequence[1]) and
            (pam_sequence[2] == 'N' or potential_pam[2] == pam_sequence[2])
        )
        
        if is_pam:
            gRNA_spacer = dna_seq[i : i + spacer_length]
            gRNAs.append({
                'spacer': str(gRNA_spacer),
                'pam': str(potential_pam),
                'strand': '+',
                'start': i,
                'end': i + spacer_length
            })

    # 2. Search for gRNAs on the reverse complement strand
    rev_comp_dna = dna_seq.reverse_complement()
    
    for i in range(len(rev_comp_dna) - spacer_length - len(pam_sequence) + 1):
        potential_pam = rev_comp_dna[i + spacer_length : i + spacer_length + len(pam_sequence)]
        is_pam = (
            (pam_sequence[0] == 'N' or potential_pam[0] == pam_sequence[0]) and
            (pam_sequence[1] == 'N' or potential_pam[1] == pam_sequence[1]) and
            (pam_sequence[2] == 'N' or potential_pam[2] == pam_sequence[2])
        )
        
        if is_pam:
            gRNA_spacer = rev_comp_dna[i : i + spacer_length]
            gRNAs.append({
                'spacer': str(gRNA_spacer),
                'pam': str(potential_pam),
                'strand': '-',
                'start': len(dna_seq) - (i + spacer_length + len(pam_sequence)),
                'end': len(dna_seq) - (i + len(pam_sequence))
            })

    return gRNAs

def score_gRNA(gRNA_info):
    """
    A very simple example of an on-target scoring function with error handling.
    """
    spacer = gRNA_info['spacer']
    score = 0
    
    try:
        gc_content = gc_fraction(spacer) * 100
        score = gc_content
    
        if spacer[0] != 'G':
            score -= 10
    except Exception as e:
        print(f"WARNING: Could not score gRNA {spacer}. Error: {e}")
        
    gRNA_info['score'] = score
    return gRNA_info


def find_and_score_gRNAs(dna_sequence, organism="human"):
    """
    Combines the finding, scoring, and ranking of gRNAs with off-target counting.
    Now also includes the flanking nucleotide information.
    """
    genome_index_path = GENOME_INDEXES.get(organism)
    if not genome_index_path:
        raise ValueError(f"No genome index found for organism: {organism}")

    candidate_gRNAs = find_gRNAs(dna_sequence)
    
    scored_gRNAs = []

    for g in candidate_gRNAs:
        # On-target scoring
        scored_grna = score_gRNA(g)
        
        # --- NEW CODE: Extract the 3' flanking nucleotide from the original sequence ---
        # The gRNA's end position is the start of the PAM.
        # The flanking nucleotide is at the position just before the gRNA start.
        if g['strand'] == '+':
            flanking_start = g['start'] - 1
            flanking_end = g['start']
            if flanking_start >= 0:
                scored_grna['flanking_nt'] = dna_sequence[flanking_start:flanking_end]
            else:
                scored_grna['flanking_nt'] = '---'
        else:
            # For the reverse complement, we need the original sequence's flanking region
            # which is at the end of the reverse complement gRNA's position in the original sequence
            original_start = len(dna_sequence) - g['end'] - 1
            original_end = len(dna_sequence) - g['end']
            if original_start >= 0:
                scored_grna['flanking_nt'] = dna_sequence[original_start:original_end]
            else:
                scored_grna['flanking_nt'] = '---'
        
        # Off-target search to get all mismatch counts
        mismatch_counts = run_bowtie_search(scored_grna['spacer'], genome_index_path)
        
        # Count off-targets for each mismatch level
        mm_counts = {'MM0': 0, 'MM1': 0, 'MM2': 0, 'MM3': 0}
        for count in mismatch_counts:
            if count <= 3:
                mm_key = f'MM{count}'
                if mm_key in mm_counts:
                    mm_counts[mm_key] += 1
        
        scored_grna.update(mm_counts)
        
        scored_gRNAs.append(scored_grna)

    # Sort by on-target score, then by off-target counts
    ranked_gRNAs = sorted(scored_gRNAs, 
                          key=lambda x: (x['score'], -x['MM0'], -x['MM1'], -x['MM2'], -x['MM3']), 
                          reverse=True)
    
    return ranked_gRNAs


# The `if __name__ == "__main__":` block is for running this script directly
# from the command line. You can uncomment it if you want to test the script
# on its own, but it's not needed for the web app.

# if __name__ == "__main__":
#     target_sequence = "GACACCGCACGGGAGATCTTGAGAGCTACATTCTCGATGTCATGCATGGGCGAGTTCGGCCATTGTGTGAGTG"
#     ranked_gRNAs = find_and_score_gRNAs(target_sequence)
    
#     print(f"Found {len(ranked_gRNAs)} potential gRNA candidates.")
#     print("\nAll ranked gRNAs:")
    
#     for i, gRNA in enumerate(ranked_gRNAs):
#         print(f"Rank {i+1}:")
#         print(f"  Spacer: {gRNA['spacer']}")
#         print(f"  PAM: {gRNA['pam']}")
#         print(f"  Strand: {gRNA['strand']}")
#         print(f"  Coordinates: {gRNA['start']}-{gRNA['end']}")
#         print(f"  Score (simple): {gRNA['score']:.2f}")
#         print("-" * 20)
