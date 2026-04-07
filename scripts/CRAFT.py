import os
import re
import platform
import shutil
import subprocess
import webbrowser
import pandas as pd
from PIL import Image
import glob
from concurrent.futures import ProcessPoolExecutor

# --- Configuration ---
SCAFFOLD_SEQUENCE = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGC"

# --- Helper function for reverse complement ---
def reverse_complement(dna_seq):
    """Returns the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATGC", "TACG")
    return dna_seq.upper().translate(complement_map)[::-1]

# --- Helper Function for Ghostscript ---
def get_ghostscript_command():
    """Finds the correct Ghostscript command for the operating system."""
    if platform.system() == "Windows":
        gs_path = None
        for path in [os.environ.get("ProgramFiles", ""), os.environ.get("ProgramFiles(x86)", "")]:
            if path and os.path.exists(path):
                gs_dir = os.path.join(path, "gs")
                if os.path.exists(gs_dir):
                    versions = [d for d in os.listdir(gs_dir) if d.startswith('gs')]
                    if versions:
                        latest_version = sorted(versions)[-1]
                        gs_exe = os.path.join(gs_dir, latest_version, "bin", "gswin64c.exe")
                        if os.path.exists(gs_exe):
                            gs_path = gs_exe
                            break
        return gs_path
    else:
        return "gs"

def generate_colored_plots(rnafold_output, plots_dir, id_suffix=""):
    """
    Parses RNAfold output and uses RNAplot to generate colored structure plots.
    This version includes robust parsing to handle cases where RNAfold might
    fail to produce a structure for a specific sequence.
    """
    print(f"⏳ Generating colored PostScript files for '{id_suffix or 'adjusted'}' sequences...")
    lines = rnafold_output.strip().split('\n')
    for i in range(0, len(lines), 2):
        if i + 1 >= len(lines):
            continue
        header = lines[i]
        seq_struct_line = lines[i+1]
        if not header.startswith('>'):
            continue
        grna_id_with_suffix = header[1:].strip()
        if ' ' not in seq_struct_line or '(' not in seq_struct_line:
            #print(f"⚠️ Warning: Skipping plot for '{grna_id_with_suffix}' because its RNAfold output line is malformed or missing a structure.")
            #print(f"   L> Offending line: \"{seq_struct_line}\"")
            continue
        try:
            sequence_and_structure, mfe = seq_struct_line.strip().rsplit(' ', 1)
        except ValueError:
            print(f"⚠️ Warning: Skipping plot for '{grna_id_with_suffix}' due to unexpected output format.")
            print(f"   L> Offending line: \"{seq_struct_line}\"")
            continue
        sequence = "".join(filter(str.isalpha, sequence_and_structure))
        structure_string = "".join(filter(lambda char: char in '().', sequence_and_structure))
        if not sequence or not structure_string:
            print(f"⚠️ Warning: Failed to parse sequence or structure for '{grna_id_with_suffix}'. Skipping plot.")
            continue
        rnaplot_input = f">{grna_id_with_suffix}\n{sequence}\n{structure_string}\n"
        try:
            command = ["RNAplot"]
            process = subprocess.run(
                command,
                input=rnaplot_input.encode('utf-8'),
                check=True,
                cwd=output_dir,
                capture_output=True,
                text=True
            )
        except subprocess.CalledProcessError as e:
            print(f"\n❌ RNAplot command failed for '{grna_id_with_suffix}'. Make sure ViennaRNA is installed correctly.")
            print(f"   L> RNAplot stderr: {e.stderr}")
            continue

# <<< NEW SECTION: Function to generate Excel Report >>>
def generate_excel_report(gRNA_data, output_dir):
    """Generates a formatted Excel report with embedded images."""
    print("✍️ Generating Excel report...")
    plots_dir = os.getcwd()
    print("plots directory for report:", plots_dir)

    report_path = os.path.join(output_dir, "gRNA_structures_report.xlsx")

    # Create a DataFrame from the results
    df = pd.DataFrame(gRNA_data)

    # Define column order and add placeholders for images
    df['Predicted Structure Image (Original)'] = ""
    df['Predicted Structure Image (Adjusted)'] = ""
    
    # Rename columns for clarity in the report
    column_map = {
        'id': 'ID', 'original_seq': 'Original gRNA Sequence (20nt)', 'display_seq': 'Adjusted gRNA Sequence (with G at 5-end if needed',
        'pam_seq': 'PAM', 'flanking_nuc': 'Flanking nt', 'strand': 'Strand',
        'gc_content': '%GC (6nt pre-PAM)', 'mfe_original': 'MFE (Original Seq) (kcal/mol)',
        'mfe_adjusted': 'MFE (Adjusted Seq) (kcal/mol)', 'score': 'MIT Specificity Score',
        'ends_in_ctt': 'CTT at -3 to -1 (CTT NGG)', 'has_ttt': 'TTT in guide (NNTTTNN)',
        'flanking_is_c': 'C downstream of PAM (NGG C)', 'g_at_pos_20': 'G at -1 of PAM (G NGG)',
        'g_at_pos_19': 'G at -2 of PAM (GN NGG)', 'c_at_pos_18': 'C at -3 of PAM (CNN NGG)',
        'g_at_pos_17': 'G at -4 of PAM (GNNN NGG)', 'flanking_is_g': 'G downstream of PAM (NGG G)',
        't_at_pam_neg_3': 'T at N position of of PAM (NGG)', 'c_or_t_at_pos_20': 'C or T at -1 of PAM (C/T NGG)',
        'g_at_pos_1_or_gg': 'G or GG at 5-end of guide (Original Seq)'
    }
    df = df.rename(columns=column_map)

    final_columns = list(column_map.values())
    # Manually insert the image columns into the desired position
    final_columns.insert(7, 'Predicted Structure Image (Original)')
    final_columns.insert(9, 'Predicted Structure Image (Adjusted)')
    df = df[final_columns]

    # Create an Excel writer object
    with pd.ExcelWriter(report_path, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name='gRNA Report', index=False)
        workbook = writer.book
        worksheet = writer.sheets['gRNA Report']

        # --- Formatting ---
        header_format = workbook.add_format({'bold': True, 'text_wrap': True, 'valign': 'top', 'fg_color': '#D7E4BC', 'border': 1})
        good_format = workbook.add_format({'bg_color': '#D4EDDA', 'font_color': '#155724'})
        bad_format = workbook.add_format({'bg_color': '#F8D7DA', 'font_color': '#721C24'})
        
        # Write header with format
        for col_num, value in enumerate(df.columns.values):
            worksheet.write(0, col_num, value, header_format)

        # Apply conditional formatting
        # Note: Column letters might change if you reorder. Adjust if necessary.
        rules = {
            'M': ('Yes', bad_format, good_format),   # CTT at -3 to -1 -> Bad
            'N': ('Yes', bad_format, good_format),   # TTT in guide -> Bad
            'O': ('Yes', good_format, bad_format),  # C downstream -> Good
            'P': ('Yes', good_format, bad_format),  # G at -1 -> Good
            'Q': ('Yes', good_format, bad_format),  # G at -2 -> Good
            'R': ('Yes', good_format, bad_format),  # C at -3 -> Good
            'S': ('Yes', good_format, bad_format),  # G at -4 -> Good
            'T': ('Yes', bad_format, good_format),   # G downstream -> Bad
            'U': ('Yes', bad_format, good_format),   # T at N -> Bad
            'V': ('Yes', bad_format, good_format),   # C or T at -1 -> Bad
            'W': ('Yes', good_format, bad_format)   # G or GG at 5' -> Good
        }
        for col, (val, fmt_yes, fmt_no) in rules.items():
            worksheet.conditional_format(f'{col}2:{col}{len(df) + 1}', {'type': 'cell', 'criteria': '==', 'value': f'"{val}"', 'format': fmt_yes})
            worksheet.conditional_format(f'{col}2:{col}{len(df) + 1}', {'type': 'cell', 'criteria': '!=', 'value': f'"{val}"', 'format': fmt_no})


        # --- Insert Images and Set Row/Column Sizes ---
        worksheet.set_column('B:C', 25)  # Sequence columns
        worksheet.set_column('H:H', 30)  # Original Image column
        worksheet.set_column('J:J', 30)  # Adjusted Image column

        # Image insertion options
        image_options = {
            'object_position': 1,
            'x_offset': 0,
            'y_offset': 0,
            'width': 600,
            'height': 150
        }

        for idx, grna in enumerate(gRNA_data):
            row_num = idx + 1 # +1 to account for header row
            worksheet.set_row(row_num, 125) # Set row height in points
            seq_name = f" gRNA {idx + 1:03}"  # Add a space for better visual separation
            worksheet.write(row_num, 0, seq_name) # Write sequence name in the first column for clarity

            # Insert original image
            img_path_orig = os.path.join(plots_dir, f"{grna['id']}_original_ss.png")
            if os.path.exists(img_path_orig):
                worksheet.insert_image(row_num, 7, img_path_orig, image_options)

            # Insert adjusted image
            img_path_adj = os.path.join(plots_dir, f"{grna['id']}_ss.png")
            if os.path.exists(img_path_adj):
                worksheet.insert_image(row_num, 9, img_path_adj, image_options)
    
    return report_path



# --- Main Script ---
def main(excel_path):
    print(f"\nProcessing file: {excel_path}")

    try:
        df = pd.read_excel(excel_path, header=None)
    except Exception as e:
        print(f"\n❌ Error reading Excel file: {e}")
        return
    
    try:
        full_dna_sequence = str(df.iloc[1, 1]).strip().upper()
    except Exception as e:
        print(f"\n❌ Error extracting exon sequence from cell B2: {e}")
        return
    
    rev_comp_full_dna = reverse_complement(full_dna_sequence)

    target_col_index = None
    score_col_index = None
    header_row_index = None
    for r_idx, row in df.iterrows():
        row_list = list(row)
        if "targetSeq" in row_list:
            target_col_index = row_list.index("targetSeq")
            header_row_index = r_idx
        if "mitSpecScore" in row_list:
            score_col_index = row_list.index("mitSpecScore")
            if header_row_index is not None:
                break
    if target_col_index is None:
        print("\n❌ Error: Column 'targetSeq' not found in the Excel file.")
        return
    if score_col_index is None:
        print("\n❌ Error: Column 'mitSpecScore' not found in the Excel file.")
        return
    print(f"✅ Found 'targetSeq' and 'mitSpecScore' columns.")
    gRNA_data = []
    for index, row_data in df.iloc[header_row_index + 1:].iterrows():
        cell_text = str(row_data[target_col_index]).strip().upper()
        score = str(row_data[score_col_index]).strip()
        if len(cell_text) >= 23:
            original_guide_seq = cell_text[:20]
            pam_seq = cell_text[20:23]
            target_with_pam = original_guide_seq + pam_seq
            if all(c in 'ACGT' for c in original_guide_seq):
                flanking_nuc = "N/A"
                strand = "N/A"
                start_index = full_dna_sequence.find(target_with_pam)
                if start_index != -1:
                    strand = "+"
                    plus_one_index = start_index + 23
                    if plus_one_index < len(full_dna_sequence):
                        flanking_nuc = full_dna_sequence[plus_one_index]
                else:
                    start_index_rev = rev_comp_full_dna.find(target_with_pam)
                    if start_index_rev != -1:
                        strand = "-"
                        plus_one_index_rev = start_index_rev + 23
                        if plus_one_index_rev < len(rev_comp_full_dna):
                            flanking_nuc = rev_comp_full_dna[plus_one_index_rev]
                adjusted_guide_seq = original_guide_seq
                display_guide_seq = original_guide_seq
                if not original_guide_seq.startswith('G'):
                    if len(original_guide_seq) > 1 and original_guide_seq[1] == 'G':
                        display_guide_seq = f"({original_guide_seq[0]}){original_guide_seq[1:]}"
                        adjusted_guide_seq = original_guide_seq[1:]
                    else:
                        adjusted_guide_seq = 'G' + original_guide_seq[1:]
                        display_guide_seq = f"G{original_guide_seq[1:]}"
                g_at_pos_1_or_gg = "Yes" if original_guide_seq.startswith('G') else "No"
                g_at_pos_20 = "Yes" if original_guide_seq[19] == 'G' else "No"
                g_at_pos_19 = "Yes" if original_guide_seq[18] == 'G' else "No"
                c_at_pos_18 = "Yes" if original_guide_seq[17] == 'C' else "No"
                g_at_pos_17 = "Yes" if original_guide_seq[16] == 'G' else "No"
                c_or_t_at_pos_20 = "Yes" if original_guide_seq[19] in ['C', 'T'] else "No"
                t_at_pam_neg_3 = "Yes" if cell_text[20] == 'T' else "No"
                ends_in_ctt = "Yes" if original_guide_seq.endswith("CTT") else "No"
                has_ttt = "Yes" if "TTT" in original_guide_seq else "No"
                flanking_is_c = "Yes" if flanking_nuc == 'C' else "No"
                flanking_is_g = "Yes" if flanking_nuc == 'G' else "No"
                seq_for_gc = original_guide_seq[-6:]
                g_count = seq_for_gc.count('G')
                c_count = seq_for_gc.count('C')
                gc_content = ((g_count + c_count) / 6) * 100 if len(seq_for_gc) == 6 else 0
                gRNA_data.append({
                    'id': f'gRNA_{index+1}', 'original_seq': original_guide_seq, 'seq': adjusted_guide_seq,
                    'display_seq': display_guide_seq, 'score': score, 'pam_seq': pam_seq,
                    'flanking_nuc': flanking_nuc, 'strand': strand, 'gc_content': f"{gc_content:.2f}%",
                    'flanking_is_c': flanking_is_c, 'flanking_is_g': flanking_is_g, 'ends_in_ctt': ends_in_ctt,
                    'has_ttt': has_ttt, 'g_at_pos_20': g_at_pos_20, 'g_at_pos_19': g_at_pos_19,
                    'c_at_pos_18': c_at_pos_18, 'g_at_pos_17': g_at_pos_17, 't_at_pam_neg_3': t_at_pam_neg_3,
                    'c_or_t_at_pos_20': c_or_t_at_pos_20, 'g_at_pos_1_or_gg': g_at_pos_1_or_gg
                })
            else:
                print(f"⚠️ Warning: Skipping invalid extracted sequence in row {index+1}: '{original_guide_seq}'")
        else:
            print(f"⚠️ Warning: Skipping entry in row {index+1} because it has fewer than 23 characters.")
    if not gRNA_data:
        print("\n❌ No valid gRNA sequences could be extracted.")
        return
    
    # Move to output directory
    os.chdir("../output")
    #Make subdirectory for this excel file
    os.makedirs(os.path.basename(excel_path), exist_ok=True)
    os.chdir(os.path.basename(excel_path))
    output_dir = os.getcwd()
    print(f"Output directory: {output_dir}")

    #Make plots directory
    os.makedirs("plots", exist_ok=True)
    os.chdir("plots/")
    print(f"\n✅ Found {len(gRNA_data)} valid gRNAs. Saving results in 'output' directory.")

    plots_dir = os.getcwd()
    print("plots directory:", plots_dir)
    print("Working directory:", os.getcwd())

    fasta_content_adjusted = ""
    for grna in gRNA_data:
        rna_seq_adjusted = grna['seq'].replace('T', 'U')
        full_rna_seq_adjusted = rna_seq_adjusted + SCAFFOLD_SEQUENCE
        fasta_content_adjusted += f">{grna['id']}\n{full_rna_seq_adjusted}\n"
    try:
        command = ["RNAfold", "-p", "-T", "37"]
        completed_process_adjusted = subprocess.run(command, input=fasta_content_adjusted.encode('utf-8'), check=True, cwd=os.getcwd(), capture_output=True)
        rnafold_output_adjusted = completed_process_adjusted.stdout.decode('utf-8')
        print("✅ RNAfold analysis complete for adjusted sequences.")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("\n❌ Error: RNAfold command failed for adjusted sequences. Make sure ViennaRNA is installed correctly and in your system's PATH.")
        return
    mfe_pattern = re.compile(r'\(\s*(-?\d+\.\d+)\s*\)')
    mfe_values_adjusted = mfe_pattern.findall(rnafold_output_adjusted)
    if len(mfe_values_adjusted) == len(gRNA_data):
        for i, grna in enumerate(gRNA_data):
            grna['mfe_adjusted'] = mfe_values_adjusted[i]
    else:
        print("⚠️ Warning: Could not parse MFE values for all adjusted sequences.")
        for grna in gRNA_data:
            grna['mfe_adjusted'] = "N/A"
    fasta_content_original = ""
    for grna in gRNA_data:
        rna_seq_original = grna['original_seq'].replace('T', 'U')
        full_rna_seq_original = rna_seq_original + SCAFFOLD_SEQUENCE
        fasta_content_original += f">{grna['id']}_original\n{full_rna_seq_original}\n"
    try:
        command = ["RNAfold", "-p", "-T", "37"]
        completed_process_original = subprocess.run(command, input=fasta_content_original.encode('utf-8'), check=True, cwd=os.getcwd(), capture_output=True)
        rnafold_output_original = completed_process_original.stdout.decode('utf-8')
        print("✅ RNAfold analysis complete for original sequences.")
    except (subprocess.CalledProcessError, FileNotFoundError):
        print("\n❌ Error: RNAfold command failed for original sequences. Make sure ViennaRNA is installed correctly and in your system's PATH.")
        return
    mfe_values_original = mfe_pattern.findall(rnafold_output_original)
    if len(mfe_values_original) == len(gRNA_data):
        for i, grna in enumerate(gRNA_data):
            grna['mfe_original'] = mfe_values_original[i]
    else:
        print("⚠️ Warning: Could not parse MFE values for all original sequences.")
        for grna in gRNA_data:
            grna['mfe_original'] = "N/A"

    generate_colored_plots(rnafold_output_adjusted, plots_dir, id_suffix="_ss")
    generate_colored_plots(rnafold_output_original, plots_dir, id_suffix="_original")
    print("⏳ Converting structure files to PNG images...")
    gs_command = get_ghostscript_command()
    if not gs_command:
        print("\n❌ Error: Could not find Ghostscript. Please install it to convert images.")
        return
    for grna in gRNA_data:
        ps_file_adjusted = os.path.join(plots_dir, f"{grna['id']}_ss.ps")
        png_file_adjusted = os.path.join(plots_dir, f"{grna['id']}_ss.png")
        if os.path.exists(ps_file_adjusted):
            try:
                cmd_adjusted = [gs_command, "-sDEVICE=pngalpha", "-r150", "-o", png_file_adjusted, "-dEPSCrop", ps_file_adjusted]
                subprocess.run(cmd_adjusted, check=True, capture_output=True)
                # Resize the PNG to 215x165 after conversion
                if os.path.exists(png_file_adjusted):
                    img = Image.open(png_file_adjusted)
                    try:
                        resample_mode = Image.Resampling.LANCZOS
                    except AttributeError:
                        resample_mode = Image.ANTIALIAS  # for older Pillow versions
                    img = img.resize((215, 165), resample_mode)
                    img.save(png_file_adjusted)
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                print(f"⚠️ Warning: Could not convert {os.path.basename(ps_file_adjusted)}. Make sure Ghostscript is installed.")
        else:
            print(f"⚠️ Warning: PostScript file not found for adjusted sequence {grna['id']}. Skipping image conversion.")
        ps_file_original = os.path.join(plots_dir, f"{grna['id']}_original_ss.ps")
        png_file_original = os.path.join(plots_dir, f"{grna['id']}_original_ss.png")
        if os.path.exists(ps_file_original):
            try:
                cmd_original = [gs_command, "-sDEVICE=pngalpha", "-r150", "-o", png_file_original, "-dEPSCrop", ps_file_original]
                subprocess.run(cmd_original, check=True, capture_output=True)
                # Resize the PNG to 215x165 after conversion
                if os.path.exists(png_file_original):
                    img = Image.open(png_file_original)
                    try:
                        resample_mode = Image.Resampling.LANCZOS
                    except AttributeError:
                        resample_mode = Image.ANTIALIAS  # for older Pillow versions
                    img = img.resize((215, 165), resample_mode)
                    img.save(png_file_original)
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                print(f"⚠️ Warning: Could not convert {os.path.basename(ps_file_original)}. Make sure Ghostscript is installed.")
        else:
            print(f"⚠️ Warning: PostScript file not found for original sequence {grna['id']}. Skipping image conversion.")

    report_path = generate_excel_report(gRNA_data, output_dir)
    print(f"\n🎉 All done! Your Excel report is ready: {report_path}")
    os.chdir(input_dir)
    # The webbrowser line is removed as it cannot open Excel files reliably.

# --- End of Function Definintions ---

#Start of Script Execution
print("--- Running CRISPR RNA Analysis and Folding Tool ---")
    
# Change to input directory and gather all Excel file paths
if __name__ == "__main__": 
    # Change to input directory
    os.chdir("../input")
    input_dir = os.getcwd()
    excel_file_paths = [file for file in os.listdir() if file.lower().endswith(('.xls', '.xlsx'))]
    if not excel_file_paths:
        print("\n❌ Error: No valid Excel files were found. Please try again.")
        exit(1)

    print("Found the following Excel files:")
    for file in excel_file_paths:
        print(f" - {file}")

    # --- Main Execution with Parallelization ---
    try:
        # Parallel execution using ProcessPoolExecutor
        with ProcessPoolExecutor() as executor:
           executor.map(main, excel_file_paths)
    except Exception as e:
        print(f"⚠️ Parallelization failed ({e}). Running sequentially...")
        for file_path in excel_file_paths:
            main(file_path)
