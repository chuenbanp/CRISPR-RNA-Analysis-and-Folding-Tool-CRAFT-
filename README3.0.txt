# CRISPR gRNA Secondary Structure Analysis Pipeline

A Python-based pipeline for analyzing CRISPR gRNA candidates, predicting secondary structures, and generating a comprehensive Excel report with fold-structure images and design metrics.

---

### System Requirements
- **OS:** macOS, Linux, or Windows
- **Python: 3.8+ (recommended 3.9 for compatibility and conda stability)
- **Dependencies:**  
  - pandas, Pillow, xlsxwriter
  - **External tools:**  
    - https://www.tbi.univie.ac.at/RNA/ (RNAfold, RNAplot)  
    - https://www.ghostscript.com/ for image conversion  
  - Additional Excel engines:
    - openpyxl (required for .xlsx files)
    - xlrd (required for .xls files)

- **Excel input:** gRNA candidate list with columns `targetSeq` and `mitSpecScore` (oringally intended for use with CRISPOR output excel files)
    -Excel files can be .xls or .xlsx:
      - .xlsx files require openpyxl (installed via requirements.txt or conda)
      - .xls files require xlrd (installed via requirements.txt or conda)

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

3. **Check output:**  
   - Images and report will be saved in a folder named `gRNA_structures_<input_filename>`  
   - Final Excel report: `gRNA_structures_report.xlsx`

---

### Notes
- Input Excel must have `targetSeq` and `mitSpecScore` columns.
- RNAfold may take time depending on the number of gRNAs.
- Ghostscript is required for converting RNAplot `.ps` files to `.png`.


---

### Installation & Setup

#### Option 1: Using conda (recommended to avoid dependency conflicts). For an introduction to conda, please see https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html
#From within the parent CRAFT Folder 
```bash
conda env create -f CRAFT_env.yml
conda activate CRAFT_env
```

This will install Python 3.9, required libraries, and external tools if available via conda.


#### Option 2: Manually install ViennaRNA and Ghostscript, pip for rest of the required packages 
- **ViennaRNA** (RNAfold, RNAplot): https://www.tbi.univie.ac.at/RNA/
- **Ghostscript**: https://www.ghostscript.com/

#From within the parent CRAFT Folder 
```bash
pip install -r requirements.txt
```
# requirements.txt includes pandas, Pillow, XlsxWriter, openpyxl, and xlrd

---

### Portability & Troubleshooting
- The pipeline uses `subprocess` calls for ViennaRNA and Ghostscript. These must be installed and accessible.
- Python dependencies are listed in `requirements.txt` and `environment.yml`.
- If Ghostscript is missing, image conversion will fail. Consider adding a fallback or skipping image embedding.
- For Windows users, Ghostscript executable is typically `gswin64c.exe`. For Mac/Linux, it's `gs`.

---

### Runtime Dependency Checks
The script will verify:
- `RNAfold` availability
- `RNAplot` availability
- `Ghostscript` availability

If any are missing, it should print instructions and exit.
