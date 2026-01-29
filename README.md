# CRISPR RNA Analysis and Folding Tool (CRAFT)

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
1. **Input:** Excel files containing gRNA sequences and MIT specificity scores.
2. **Processing:**
   - Extracts gRNA sequences and PAM sites
   - Adds alternative sequences for U3/U6 polIII-type RNA polymerase transcription efficiency (forces leading G if needed)
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

### Installation & Setup
#### Cloning the Repository
   Download the pipeline code by cloning this GitHub repository.
   In Terminal:
   ```bash/terminal
   # Navigate to the directory where you want the repo
   cd /path/to/your/projects

   # Clone the repository
   git clone https://github.com/chuenbanp/gRNA-Design-and-Structure-Prediction.git

   # Enter the cloned folder
   cd gRNA-Design-and-Structure-Prediction/
   ```

#### Install Dependencies
   #### Option 1: Using conda (recommended to avoid dependency conflicts). 
   From within ./gRNA-Design-and-Structure-Prediction
   ```bash/terminal
   conda env create -f CRAFT_env.yml
   conda activate CRAFT_env
   ```
   This will install Python 3.9, required libraries, and external tools all available via conda.

   For an introduction to conda, please see https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html

   #### Option 2: Manually install ViennaRNA and Ghostscript, pip for rest (collected in requirements.txt)
   - **ViennaRNA** (RNAfold, RNAplot): https://www.tbi.univie.ac.at/RNA/
   - **Ghostscript**: https://www.ghostscript.com/
   From within the parent CRAFT Folder:
   ```bash/terminal
   pip install -r requirements.txt
   ```
   requirements.txt includes pandas, Pillow, XlsxWriter, openpyxl, and xlrd

---

### Running the Pipeline
   1. **Install Dependencies, see above**

   2. *Move CRISPOR/gRNA excel files into ./input directory

   3. **Run the script:**  
      Navigate to ./scripts directory
      ```bash/terminal
      cd scripts/
      "$CONDA_PREFIX/bin/python" -u CRAFT.py
      ```

   4. **Check ./output directory:**  
      - Images and report will be saved in a folders named for each input excel file
      - Final Excel report: `gRNA_structures_report.xlsx` generated in each folder per input file

---

### Notes
- Input Excel must have `targetSeq` and `mitSpecScore` columns.
- RNAfold/.png conversion may take time depending on the number of gRNAs.
- Ghostscript is required for converting RNAplot `.ps` files to `.png`.

---

### Portability & Troubleshooting
- The pipeline uses `subprocess` calls for ViennaRNA and Ghostscript. These must be installed and accessible.
- Python dependencies are listed in `requirements.txt` and `environment.yml`.
- Conda installation of required packages should be portable across Mac, Windows, and Linux.

---
