"""
Author: [Your Name]
Infrastructure: [Your Lab/Institution]
Purpose: This script processes mass spectrometry (MS) data and converts it into the MSP (Mass Spectrometry Power) format.

Overview:
---------
This code is designed to convert mass spectrometry data from a list of DataFrames into the MSP (Mass Spectrometry Power) format, typically used for storing and sharing data about metabolites and their MS/MS spectra. The script extracts key compound information, along with their associated fragment peaks, and generates an MSP file that can be used for further analysis.

Steps:
------
1. **Input**:
   - The script accepts a list of DataFrames (`dfs`), each containing mass spectrometry data with the following relevant columns:
     - `Compound`: The name of the compound.
     - `Peak Label`: A label indicating whether the entry is a main peak or a fragment peak (fragment peaks contain 'F' in the label).
     - `m/z (Expected)`: The expected m/z value of the peak.
     - `Formula`: The chemical formula of the compound.
     - `RT`: The retention time for the compound.
     - `Charge`: The ion mode or charge state of the compound.
     - `Family`: The compound class (e.g., lipid family).
     - `Adduct`: The ion adduct used in mass spectrometry.

2. **Processing**:
   - The script filters out rows where the `Peak Label` contains 'F' (indicating fragment peaks), and processes only the "main" entries.
   - For each main entry, it extracts relevant information (e.g., compound name, precursor m/z, formula, retention time, etc.).
   - The script checks if the retention time is available and non-zero, formatting it accordingly.
   - For each compound, the script then finds its associated fragment peaks (rows where the `Peak Label` contains 'F') and collects their m/z values.

3. **MSP Entry Construction**:
   - For each compound, a properly formatted MSP entry is created, including the following details:
     - `NAME`: The name of the compound.
     - `PRECURSORMZ`: The precursor m/z value.
     - `PRECURSORTYPE`: The adduct type.
     - `FORMULA`: The molecular formula of the compound.
     - `RETENTIONTIME`: The formatted retention time.
     - `IONMODE`: The ion mode (charge state).
     - `COMPOUNDCLASS`: The compound class.
     - `Comment`: A comment indicating the theoretical MS2 spectrum was created from Orbitrap Lipidomics data.
     - `Num Peaks`: The number of fragment peaks associated with the compound.
     - The fragment peaks themselves, listed as m/z values with a placeholder intensity value of `999`.

4. **Output**:
   - The script writes all the generated MSP entries to a file named `HomeDB_MS2.msp`, which will be saved in the current directory.

5. **Logging**:
   - A success message is printed once the data is successfully written to the MSP file.

Input Requirements:
-------------------
- A list of DataFrames (`dfs`), where each DataFrame contains columns for compound data and peak information as described above.

Output:
-------
- **MSP File**: An MSP-formatted text file (`HomeDB_MS2.msp`) containing the processed compound and fragment peak information.
- **Log Message**: A success message indicating the completion of the MSP file generation.

Code Modules and Libraries:
---------------------------
- pandas: For handling DataFrame operations and extracting relevant compound information.
- os: For file and directory management (optional depending on file saving logic).
- re: For filtering and processing data based on regular expressions.
- logging: For tracking and logging information about the data processing (if added).

Usage:
------
To use this script, provide a list of DataFrames (`dfs`) that contain the necessary mass spectrometry data. The script will process the data, convert it into MSP format, and output the results to a file named `HomeDB_MS2.msp`. After running the script, you will see a success message confirming the output.

Note:
- This script assumes that the compound names in each DataFrame are unique, and the fragment peaks are correctly labeled with 'F' in the `Peak Label` column.
"""

import pandas as pd

dfs = ['liste of dfs'] # put your dfs for convert them to msp

output_lines = []

for df in dfs:
    main_entries = df[~df['Peak Label'].str.contains('F', na=False)]

    for _, main_row in main_entries.iterrows():
        compound_name = main_row['Compound']
        precursor_mz = main_row['m/z (Expected)']
        formula = main_row['Formula']
        retention_time = main_row['RT']
        ion_mode = main_row['Charge']
        compound_class = main_row['Family']
        adduct = main_row['Adduct']

        retention_time_str = f"RETENTIONTIME: {retention_time}" if pd.notna(retention_time) and retention_time != 0 else "RETENTIONTIME: "

        num_peaks_rows = df[(df['Compound'] == compound_name) & (df['Peak Label'].str.contains('F', na=False))]

        num_peaks = []
        for _, peak_row in num_peaks_rows.iterrows():
            num_peaks.append(f"{peak_row['m/z (Expected)']:.5f}\t999")
        num_peaks_str = "\n".join(num_peaks)

        entry = f"""NAME: {compound_name}
PRECURSORMZ: {precursor_mz:.5f}
PRECURSORTYPE: {adduct}
SMILES: 
INCHIKEY: 
FORMULA: {formula}
{retention_time_str}
CCS: 
IONMODE: {ion_mode}
COMPOUNDCLASS: {compound_class}
Comment: theoretical MS2 created from the information of Orbitrap Lipidomics.
Num Peaks: {len(num_peaks)}
{num_peaks_str}
"""
        output_lines.append(entry)

with open("HomeDB_MS2.msp", "w") as f:
    f.write("\n".join(output_lines))

print("Data successfully written to HomeDB_MS2.msp.")