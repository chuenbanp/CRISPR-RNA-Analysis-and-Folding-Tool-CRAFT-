# app.py
import os
import re
import platform
import subprocess
import uuid
import shutil
from flask import Flask, render_template, request, send_from_directory

# --- Import your new gRNA design script ---
from gRNA_design import find_and_score_gRNAs

# --- Flask App Setup ---
app = Flask(__name__)
RESULT_FOLDER = 'results'
app.config['RESULT_FOLDER'] = RESULT_FOLDER
os.makedirs(RESULT_FOLDER, exist_ok=True)


# --- Configuration & Helpers ---
SCAFFOLD_SEQUENCE = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC"

def reverse_complement(dna_seq):
    complement_map = str.maketrans("ATGC", "TACG")
    return dna_seq.upper().translate(complement_map)[::-1]

def get_ghostscript_command():
    if platform.system() == "Windows":
        for path in [os.environ.get("ProgramFiles", ""), os.environ.get("ProgramFiles(x86)", "")]:
            if path and os.path.exists(os.path.join(path, "gs")):
                gs_dir = os.path.join(path, "gs")
                versions = sorted([d for d in os.listdir(gs_dir) if d.startswith('gs')], reverse=True)
                if versions:
                    gs_exe = os.path.join(gs_dir, versions[0], "bin", "gswin64c.exe")
                    if os.path.exists(gs_exe): return gs_exe
        return None
    else:
        return shutil.which("gs")

# --- Refactored processing function to work with gRNA_data directly ---
def analyze_grna_structures(gRNA_data, output_dir, base_name):
    """
    Takes a list of gRNA dictionaries (from the design tool) and runs
    RNAfold and Ghostscript on them.
    """
    if not gRNA_data:
        return None, "No valid gRNA sequences provided."

    # Process each gRNA dictionary to add derived parameters
    for i, grna in enumerate(gRNA_data):
        # We need to map the keys from the design tool to the keys
        # used in the original `process_files` function.
        original_guide_seq = grna['spacer']
        pam_seq = grna['pam']
        
        # Original logic for adjusting the sequence
        adjusted_guide_seq = original_guide_seq
        display_guide_seq = original_guide_seq

        if not original_guide_seq.startswith('G'):
            if len(original_guide_seq) > 1 and original_guide_seq[1] == 'G':
                display_guide_seq = f"(<span style='color:red;'>{original_guide_seq[0]}</span>){original_guide_seq[1:]}"
                adjusted_guide_seq = original_guide_seq[1:]
            else:
                adjusted_guide_seq = 'G' + original_guide_seq[1:]
                display_guide_seq = f"<span style='color:red;'>G</span>{original_guide_seq[1:]}"

        seq_for_gc = original_guide_seq[-6:]
        g_count, c_count = seq_for_gc.count('G'), seq_for_gc.count('C')

        # Add all the new keys to the dictionary for consistency
        grna['id'] = f'{base_name}_g{i+1}' # Set a unique ID
        grna['original_seq'] = original_guide_seq
        grna['adjusted_seq'] = adjusted_guide_seq
        grna['display_seq'] = display_guide_seq
        grna['pam_seq'] = pam_seq
        grna['gc_content'] = ((g_count + c_count) / 6) * 100
        grna['flanking_is_c'] = "N/A" # These parameters were from the full gDNA sequence, which we don't have anymore
        grna['flanking_is_g'] = "N/A"
        grna['ends_in_ctt'] = "Yes" if original_guide_seq.endswith("CTT") else "No"
        grna['has_ttt'] = "Yes" if "TTT" in original_guide_seq else "No"
        grna['g_at_pos_20'] = "Yes" if original_guide_seq[19] == 'G' else "No"
        grna['g_at_pos_19'] = "Yes" if original_guide_seq[18] == 'G' else "No"
        grna['c_at_pos_18'] = "Yes" if original_guide_seq[17] == 'C' else "No"
        grna['g_at_pos_17'] = "Yes" if original_guide_seq[16] == 'G' else "No"
        grna['t_at_pam_neg_3'] = "Yes" if pam_seq[0] == 'T' else "No"
        grna['c_or_t_at_pos_20'] = "Yes" if original_guide_seq[19] in ['C', 'T'] else "No"
        grna['g_at_pos_1_or_gg'] = "Yes" if adjusted_guide_seq.startswith('G') else "No"
       
        
    print(f"Running RNAfold for {len(gRNA_data)} sequences...")
    command = ["RNAfold", "-p", "--noLP", "-T", "37"]
    mfe_pattern = re.compile(r'\(\s*(-?\d+\.\d+)\s*\)')
    
    # Run RNAfold for the adjusted sequence
    fasta_content_adjusted = "".join([f">{g['id']}\n{g['adjusted_seq'].replace('T', 'U') + SCAFFOLD_SEQUENCE}\n" for g in gRNA_data])
    try:
        p_adj = subprocess.run(command, input=fasta_content_adjusted.encode('utf-8'), check=True, cwd=output_dir, capture_output=True)
        mfe_values_adjusted = mfe_pattern.findall(p_adj.stdout.decode('utf-8'))
        if len(mfe_values_adjusted) == len(gRNA_data):
            for i, grna in enumerate(gRNA_data): grna['mfe_adjusted'] = mfe_values_adjusted[i]
        else:
            for grna in gRNA_data: grna['mfe_adjusted'] = "N/A"
    except subprocess.CalledProcessError as e: return None, f"RNAfold failed for adjusted sequences. Error: {e.stderr.decode()}"
    
    # Run RNAfold for the original sequence
    fasta_content_original = "".join([f">{g['id']}_original\n{g['original_seq'].replace('T', 'U') + SCAFFOLD_SEQUENCE}\n" for g in gRNA_data])
    try:
        p_orig = subprocess.run(command, input=fasta_content_original.encode('utf-8'), check=True, cwd=output_dir, capture_output=True)
        mfe_values_original = mfe_pattern.findall(p_orig.stdout.decode('utf-8'))
        if len(mfe_values_original) == len(gRNA_data):
            for i, grna in enumerate(gRNA_data): grna['mfe_original'] = mfe_values_original[i]
        else:
            for grna in gRNA_data: grna['mfe_original'] = "N/A"
    except subprocess.CalledProcessError as e: return None, f"RNAfold failed for original sequences. Error: {e.stderr.decode()}"
    
    print("Converting PostScript files to PNG with Ghostscript...")
    gs_command = get_ghostscript_command()
    if not gs_command: return None, "Ghostscript not found. Please install it."
    
    for grna in gRNA_data:
        for suffix in ["", "_original"]:
            ps_file = os.path.join(output_dir, f"{grna['id']}{suffix}_ss.ps")
            png_file = os.path.join(output_dir, f"{grna['id']}{suffix}_ss.png") # Corrected file name
            if os.path.exists(ps_file):
                try:
                    cmd = [gs_command, "-sDEVICE=pngalpha", "-r150", "-o", png_file, "-dEPSCrop", ps_file]
                    subprocess.run(cmd, check=True, capture_output=True)
                except subprocess.CalledProcessError as e: return None, f"Ghostscript failed. Error: {e.stderr.decode()}"
    
    print("Analysis complete.")
    return gRNA_data, None


# --- Flask Routes ---
@app.route('/')
def index():
    return render_template('index.html')

@app.route('/process', methods=['POST'])
def process():
    dna_sequence = request.form.get('dna_sequence', '').strip().upper()
    organism = request.form.get('organism', 'human') # Get the user's choice
    base_name = request.form.get('base_name', 'gRNA').strip()
    
    if not dna_sequence:
        return render_template('results.html', error="Input error: Please provide a DNA sequence.")
    if not base_name:
        base_name = 'gRNA'

    session_id = str(uuid.uuid4())
    session_dir = os.path.join(app.config['RESULT_FOLDER'], session_id)
    os.makedirs(session_dir, exist_ok=True)

    print(f"\n--- New Analysis Started [Session: {session_id}] ---")
    print("Running gRNA design tool...")

    try:
        # Pass the organism variable directly to the function.
        # This function will handle accessing GENOME_INDEXES internally.
        ranked_gRNAs = find_and_score_gRNAs(dna_sequence, organism=organism)
        
        gRNA_data, error = analyze_grna_structures(ranked_gRNAs, session_dir, base_name)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return render_template('results.html', error=f"An unexpected error occurred: {e}")
    
    if error:
        print(f"An error occurred during analysis: {error}")
        return render_template('results.html', error=error)
    
    print("Rendering final results template.")
    return render_template('results.html', gRNA_data=gRNA_data, session_id=session_id)

@app.route('/results/<session_id>/<filename>')
def serve_result_file(session_id, filename):
    """This route serves the generated image files."""
    session_dir = os.path.join(app.config['RESULT_FOLDER'], session_id)
    return send_from_directory(session_dir, filename)

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
