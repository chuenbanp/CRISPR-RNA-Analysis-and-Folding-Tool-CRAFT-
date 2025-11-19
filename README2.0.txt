# CRISPR gRNA Secondary Structure Analysis Pipeline

A Python-based pipeline for analyzing CRISPR gRNA candidates, predicting secondary structures, and generating a comprehensive Excel report with embedded images and key design metrics.

---

### System Requirements
- **OS:** macOS, Linux, or Windows
- **Python:** 3.8+
- **Dependencies:**  
  - pandas, Pillow, xlsxwriter
  - **External tools:**  
    - https://www.tbi.univie.ac.at/RNA/ (RNAfold, RNAplot)  
    - https://www.ghostscript.com/ for image conversion  
- **Excel input:** gRNA candidate list with columns `targetSeq` and `mitSpecScore` (oringally intended for use with CRISPOR output excel files)

---

### Pipeline Overview
1. **Input:** Excel file containing gRNA sequences and MIT specificity scores.
2. **Processing:**  
   - Extracts gRNA sequences and PAM sites  
   - Adjusts sequences for U3/U6 polIII-type RNA polymerase transcription efficiency (adds leading G if needed)  
   - Predicts secondary structures using RNAfold  
   - Converts RNAplot PostScript outputs to PNG images  
3. **Output:**  
   - Folder with structure images  
   - Excel report summarizing:
     - Original and adjusted sequences  
     - Strand orientation and flanking nucleotides  
     - GC content, MFE values  
     - Design flags (e.g., TTT motifs, PAM context)  
     - Embedded structure images  

---

### Running the Pipeline
1. **Install dependencies:**  
   ```bash
   pip install pandas Pillow xlsxwriter
   ```
   Ensure ViennaRNA and Ghostscript are installed and available in your PATH.

2. **Run the script:**  
   ```bash
   python gRNA_Analysis_Pipeline.py
   ```
   When prompted, drag and drop your Excel file and press Enter.

3. **Check output:**  
   - Images and report will be saved in a folder named `gRNA_structures_<input_filename>`  
   - Final Excel report: `gRNA_structures_report.xlsx`

---

### Notes
- Input Excel must have `targetSeq` and `mitSpecScore` columns.
- RNAfold may take time depending on the number of gRNAs.
- Ghostscript is required for converting RNAplot `.ps` files to `.png`.
