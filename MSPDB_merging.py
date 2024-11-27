"""
Author: Ramiz KHALED <ramiz.khaled@inserm.fr>
Infrastructure: MetaToul, plateau Lipidomique

### Purpose:
This processe fusion several MSP files, it updates data in a MSP file (`PRECURSORMZ`, `RETENTIONTIME`, `PRECURSORTYPE`, and `IONMODE`) based on provided table data MSP data bases. 
It checks if a compound name exists in an initial MSP file, updates the existing entry with information from a second MSP file, 
or creates a new entry if the compound name is not found.

### Functionality:
After reading and pre-rocessing and merging of the MSP files of interest :
1. **Reading and Parsing the MSP File**:
   - The script reads an MSP file, which contains mass spectrometry data entries.
   - Each entry in the MSP file begins with a 'NAME' field, followed by key-value pairs, such as `PRECURSORMZ`, `PRECURSORTYPE`, etc.
   - The data is normalized and split into individual compound entries for processing.
   
2. **Updating MSP Entries Based on Table Data**:
   - The table provided is expected to contain columns: `NAME`, `PRECURSORMZ`, `RETENTIONTIME`, `PRECURSORTYPE`, and `IONMODE`.
   - The script updates the corresponding fields in the MSP data based on matching compound names from the table.
   - If a compound name from the table does not exist in the MSP file, it adds a new entry with data from the table, leaving fields that are not available empty.

3. **Extracting Data**:
   - The script extracts key information for each entry such as:
     - `NAME`: The compound name.
     - `PRECURSORMZ`: The precursor mass-to-charge ratio.
     - `PRECURSORTYPE`: The type of precursor ion.
     - `RETENTIONTIME`: The retention time in the chromatographic run.
     - `IONMODE`: Ionization mode (e.g., Positive or Negative).
     - Additional fields like `SMILES`, `INCHIKEY`, `FORMULA`, `CCS`, `COMPOUNDCLASS`, `Comment`, and `Num Peaks`.
   - If any of these fields are missing in the MSP file, they are either skipped or left empty when adding or updating data.

4. **Rebuilding the MSP File**:
   - After updating or adding new entries, the script rebuilds the MSP file.
   - The new content is written to a new file (`updated_msp_file.msp`), preserving the format of the original MSP file.

5. **Output**:
   - The modified MSP file is saved as `updated_msp_file.msp` in the same directory as the script.
   - The output retains the original MSP format while incorporating the updated or newly added data.

### Requirements:
- The input MSP file must be properly formatted with fields such as `NAME`, `PRECURSORMZ`, `PRECURSORTYPE`, etc.
- The table provided should be in a Pandas DataFrame with columns `NAME`, `PRECURSORMZ`, `RETENTIONTIME`, `PRECURSORTYPE`, and `IONMODE`.
- The compound names in the MSP file and the table should match exactly or after normalization (e.g., case-insensitive).

### Example Usage:
```python
# Sample table data in a DataFrame
table_data = {
    'NAME': ['CerP 19:1;2O/26:1', 'Cholesterol'],
    'PRECURSORMZ': [770.64277, 409.34409],
    'RETENTIONTIME': [13.11, 9.05],
    'PRECURSORTYPE': ['[M+H]+', '[M+Na]+'],
    'IONMODE': ['Positive', 'Positive']
}
table_df = pd.DataFrame(table_data)

# Path to the existing MSP file (homDB_path)
homDB_path = 'homDB.msp'  # Input MSP file path

# Call the function to update the MSP file with table data
modify_msp_file(homDB_path, table_df)

"""
import pandas as pd
import numpy as np
import re
import os

homeDB_path = "PATH/FOR/YOUR/MSP1.msp"
POS_msp = "PATH/FOR/YOUR/MSP2.msp"
NEG_msp = "PATH/FOR/YOUR/MSP3.msp"

###Data reading

#Read MSP2

with open(homeDB_path, "r") as file:
    homDB_text = file.read()

normalized_text = re.sub(r'\s+', ' ', homDB_text)

entries = normalized_text.split("NAME: ")[1:]  

matches = []

for entry in entries:
    name_match = re.search(r"^(.+?)\s+PRECURSORMZ:", entry)
    if name_match:
        compound_name = name_match.group(1)
    else:
        compound_name = "N/A"

    mz_match = re.search(r"PRECURSORMZ:\s*(\S+)", entry)
    if mz_match:
        precursor_mz = mz_match.group(1)
    else:
        precursor_mz = "N/A"

    adduct_match = re.search(r"PRECURSORTYPE:\s*(\S+)", entry)
    if adduct_match:
        adduct = adduct_match.group(1)
    else:
        adduct = "N/A"

    retention_match = re.search(r"RETENTIONTIME:\s*(\S+)", entry)
    if retention_match:
        retention_time = retention_match.group(1)
    else:
        retention_time = "N/A"

    ionmode_match = re.search(r"IONMODE:\s*(\S+)", entry)
    if ionmode_match:
        ion_mode = ionmode_match.group(1)
    else:
        ion_mode = "N/A"

    matches.append((compound_name, precursor_mz, adduct, retention_time, ion_mode))

if matches:
    homDB_data = pd.DataFrame(matches, columns=["NAME", "PRECURSORMZ", "PRECURSORTYPE",  "RETENTIONTIME", "IONMODE"])
    print(homDB_data.head())
else:
    print("No matches found.")

#Read MSP2

with open(POS_msp, "r") as file:
    POS_txt = file.read()

normalized_text = re.sub(r'\s+', ' ', POS_txt)

entries = normalized_text.split("NAME: ")[1:]  

matches = []

for entry in entries:
    name_match = re.search(r"^(.+?)\s+PRECURSORMZ:", entry)
    if name_match:
        compound_name = name_match.group(1)
    else:
        compound_name = "N/A"

    mz_match = re.search(r"PRECURSORMZ:\s*(\S+)", entry)
    if mz_match:
        precursor_mz = mz_match.group(1)
    else:
        precursor_mz = "N/A"

    adduct_match = re.search(r"PRECURSORTYPE:\s*(\S+)", entry)
    if adduct_match:
        adduct = adduct_match.group(1)
    else:
        adduct = "N/A"

    retention_match = re.search(r"RETENTIONTIME:\s*(\S+)", entry)
    if retention_match:
        retention_time = retention_match.group(1)
    else:
        retention_time = "N/A"

    ionmode_match = re.search(r"IONMODE:\s*(\S+)", entry)
    if ionmode_match:
        ion_mode = ionmode_match.group(1)
    else:
        ion_mode = "N/A"

    matches.append((compound_name, precursor_mz, adduct, retention_time, ion_mode))

if matches:
    POS_data = pd.DataFrame(matches, columns=["NAME", "PRECURSORMZ", "PRECURSORTYPE", "RETENTIONTIME", "IONMODE"])
    print(POS_data.head())
else:
    print("No matches found.")

#Read MSP3

with open(NEG_msp, "r") as file:
    NEG_txt = file.read()

normalized_text = re.sub(r'\s+', ' ', NEG_txt)

entries = normalized_text.split("NAME: ")[1:]  

matches = []

for entry in entries:
    name_match = re.search(r"^(.+?)\s+PRECURSORMZ:", entry)
    if name_match:
        compound_name = name_match.group(1)
    else:
        compound_name = "N/A"

    mz_match = re.search(r"PRECURSORMZ:\s*(\S+)", entry)
    if mz_match:
        precursor_mz = mz_match.group(1)
    else:
        precursor_mz = "N/A"

    adduct_match = re.search(r"PRECURSORTYPE:\s*(\S+)", entry)
    if adduct_match:
        adduct = adduct_match.group(1)
    else:
        adduct = "N/A"

    retention_match = re.search(r"RETENTIONTIME:\s*(\S+)", entry)
    if retention_match:
        retention_time = retention_match.group(1)
    else:
        retention_time = "N/A"

    ionmode_match = re.search(r"IONMODE:\s*(\S+)", entry)
    if ionmode_match:
        ion_mode = ionmode_match.group(1)
    else:
        ion_mode = "N/A"

    matches.append((compound_name, precursor_mz, adduct, retention_time, ion_mode))

if matches:
    NEG_data = pd.DataFrame(matches, columns=["NAME", "PRECURSORMZ", "PRECURSORTYPE", "RETENTIONTIME", "IONMODE"])
    print(NEG_data.head())
else:
    print("No matches found.")

### Pre-Processing

# For MSP1

homDB_data['NAME'] = homDB_data['NAME'].str.replace(r'^Carnitine','Car', regex = True)
homDB_data['PRECURSORTYPE'] = homDB_data['PRECURSORTYPE'].where(homDB_data['PRECURSORTYPE'].str.startswith('[M'), np.nan)
homDB_data['RETENTIONTIME'] = pd.to_numeric(homDB_data['RETENTIONTIME'], errors='coerce')
homDB_data['NAME'] = homDB_data['NAME'].str.lower()

# For MSP2

POS_data['NAME'] = POS_data['NAME'].str.replace(r'^CAR','Car', regex = True)
POS_data['PRECURSORTYPE'] = POS_data['PRECURSORTYPE'].where(POS_data['PRECURSORTYPE'].str.startswith('[M'), np.nan)
POS_data['RETENTIONTIME'] = pd.to_numeric(POS_data['RETENTIONTIME'], errors='coerce')
POS_data['NAME'] = POS_data['NAME'].str.lower()

# For MSP3

NEG_data['NAME'] = NEG_data['NAME'].str.replace(r'^CAR','Car', regex = True)
NEG_data['PRECURSORTYPE'] = NEG_data['PRECURSORTYPE'].where(NEG_data['PRECURSORTYPE'].str.startswith('[M'), np.nan)
NEG_data['RETENTIONTIME'] = pd.to_numeric(NEG_data['RETENTIONTIME'], errors='coerce')
NEG_data['NAME'] = NEG_data['NAME'].str.lower()

### Merging

# For MSP1 + MSP2
#Priority for the MSP1

merged_data = pd.merge(
    POS_data, 
    homDB_data, 
    on='NAME', 
    how='outer', 
    suffixes=('_POS', '_HOM')
)
merged_data

merged_data['PRECURSORMZ'] = merged_data['PRECURSORMZ_HOM'].where(merged_data['PRECURSORMZ_HOM'].notna(), merged_data['PRECURSORMZ_POS'])

merged_data['PRECURSORTYPE'] = merged_data['PRECURSORTYPE_HOM']
merged_data['RETENTIONTIME'] = merged_data['RETENTIONTIME_HOM']

merged_data['IONMODE'] = merged_data['IONMODE_HOM'].where(merged_data['IONMODE_HOM'].notna(), merged_data['IONMODE_POS'])

Home_POS_data = merged_data[['NAME', 'PRECURSORMZ', 'RETENTIONTIME', 'PRECURSORTYPE','IONMODE']]

# For (MSP1 + MSP2) + MSP3
#Priority for the (MSP1 + MSP2)

merged_data = pd.merge(
    NEG_data, 
    Home_POS_data, 
    on='NAME', 
    how='outer', 
    suffixes=('_NEG', '_HOM')
)
merged_data

merged_data['PRECURSORMZ'] = merged_data['PRECURSORMZ_HOM'].where(merged_data['PRECURSORMZ_HOM'].notna(), merged_data['PRECURSORMZ_NEG'])

merged_data['PRECURSORTYPE'] = merged_data['PRECURSORTYPE_HOM']
merged_data['RETENTIONTIME'] = merged_data['RETENTIONTIME_HOM']

merged_data['IONMODE'] = merged_data['IONMODE_HOM'].where(merged_data['IONMODE_HOM'].notna(), merged_data['IONMODE_NEG'])

total_data = merged_data[['NAME', 'PRECURSORMZ', 'RETENTIONTIME', 'PRECURSORTYPE','IONMODE']]

### New MSP generating

def modify_txt_file(homDB_path, table_df):

    with open(homDB_path, "r") as file:
        homDB_text = file.read()

    normalized_text = re.sub(r'\s+', ' ', homDB_text)
    entries = normalized_text.split("NAME: ")[1:]  

    txt_data = {}

    current_name = None
    for entry in entries:
        lines = entry.split("\n")

        for line in lines:
            if line.startswith('NAME:'):
                current_name = line.strip().split("NAME: ")[1]
                txt_data[current_name] = {
                    'PRECURSORMZ': None,
                    'PRECURSORTYPE': None,
                    'RETENTIONTIME': None,
                    'SMILES': None,
                    'INCHIKEY': None,
                    'FORMULA': None,
                    'CCS': None,
                    'IONMODE': None,
                    'COMPOUNDCLASS': None,
                    'Comment': None,
                    'Num Peaks': 0,
                    'PEAKS': []
                }

            elif line.startswith('PRECURSORMZ:'):
                if current_name:
                    txt_data[current_name]['PRECURSORMZ'] = line.strip().split("PRECURSORMZ: ")[1]

            elif line.startswith('PRECURSORTYPE:'):
                if current_name:
                    txt_data[current_name]['PRECURSORTYPE'] = line.strip().split("PRECURSORTYPE: ")[1]

            elif line.startswith('RETENTIONTIME:'):
                if current_name:
                    txt_data[current_name]['RETENTIONTIME'] = line.strip().split("RETENTIONTIME: ")[1]

            elif line.startswith('SMILES:'):
                if current_name:
                    smi = line.strip().split("SMILES: ")[1] if len(line.strip().split("SMILES: ")) > 1 else ""
                    txt_data[current_name]['SMILES'] = smi

            elif line.startswith('INCHIKEY:'):
                if current_name:
                    inchi = line.strip().split("INCHIKEY: ")[1] if len(line.strip().split("INCHIKEY: ")) > 1 else ""
                    txt_data[current_name]['INCHIKEY'] = inchi

            elif line.startswith('FORMULA:'):
                if current_name:
                    formula = line.strip().split("FORMULA: ")[1] if len(line.strip().split("FORMULA: ")) > 1 else ""
                    txt_data[current_name]['FORMULA'] = formula

            elif line.startswith('CCS:'):
                if current_name:
                    ccs = line.strip().split("CCS: ")[1] if len(line.strip().split("CCS: ")) > 1 else ""
                    txt_data[current_name]['CCS'] = ccs

            elif line.startswith('IONMODE:'):
                if current_name:
                    ion_mode = line.strip().split("IONMODE: ")[1] if len(line.strip().split("IONMODE: ")) > 1 else ""
                    txt_data[current_name]['IONMODE'] = ion_mode

            elif line.startswith('COMPOUNDCLASS:'):
                if current_name:
                    compound_class = line.strip().split("COMPOUNDCLASS: ")[1] if len(line.strip().split("COMPOUNDCLASS: ")) > 1 else ""
                    txt_data[current_name]['COMPOUNDCLASS'] = compound_class

            elif line.startswith('Comment:'):
                if current_name:
                    comment = line.strip().split("Comment: ")[1] if len(line.strip().split("Comment: ")) > 1 else ""
                    txt_data[current_name]['Comment'] = comment

            elif line.startswith('Num Peaks:'):
                if current_name:
                    num_peaks = int(line.strip().split("Num Peaks: ")[1]) if len(line.strip().split("Num Peaks: ")) > 1 else 0
                    txt_data[current_name]['Num Peaks'] = num_peaks

            if current_name and line.strip() and line.strip().count("\t") == 1:
                txt_data[current_name]['PEAKS'].append(line.strip())

    for idx, row in table_df.iterrows():
        name_in_table = row['NAME']
        if name_in_table in txt_data:
            txt_data[name_in_table]['PRECURSORMZ'] = row['PRECURSORMZ']
            txt_data[name_in_table]['RETENTIONTIME'] = row['RETENTIONTIME']
            txt_data[name_in_table]['PRECURSORTYPE'] = row['PRECURSORTYPE']
            txt_data[name_in_table]['IONMODE'] = row['IONMODE']
        else:
            txt_data[name_in_table] = {
                'PRECURSORMZ': row['PRECURSORMZ'],
                'RETENTIONTIME': row['RETENTIONTIME'],
                'PRECURSORTYPE': row['PRECURSORTYPE'],
                'IONMODE': row['IONMODE'],
                'SMILES': "",
                'INCHIKEY': "",
                'FORMULA': "",
                'CCS': "",
                'COMPOUNDCLASS': "",
                'Comment': "",
                'Num Peaks': 0,
                'PEAKS': []
            }

    updated_txt_lines = []
    for name, data in txt_data.items():
        updated_txt_lines.append(f"NAME: {name}")
        updated_txt_lines.append(f"PRECURSORMZ: {data['PRECURSORMZ']}")
        updated_txt_lines.append(f"PRECURSORTYPE: {data['PRECURSORTYPE']}")
        updated_txt_lines.append(f"SMILES: {data['SMILES']}")
        updated_txt_lines.append(f"INCHIKEY: {data['INCHIKEY']}")
        updated_txt_lines.append(f"FORMULA: {data['FORMULA']}")
        updated_txt_lines.append(f"RETENTIONTIME: {data['RETENTIONTIME']}")
        updated_txt_lines.append(f"CCS: {data['CCS']}")
        updated_txt_lines.append(f"IONMODE: {data['IONMODE']}")
        updated_txt_lines.append(f"COMPOUNDCLASS: {data['COMPOUNDCLASS']}")
        updated_txt_lines.append(f"Comment: {data['Comment']}")
        updated_txt_lines.append(f"Num Peaks: {data['Num Peaks']}")

        for peak in data['PEAKS']:
            updated_txt_lines.append(peak)

        updated_txt_lines.append("") 

    current_directory = os.getcwd()
    output_path = os.path.join(current_directory, "updated_homeDB.msp")

    with open(output_path, 'w') as file:
        file.writelines("\n".join(updated_txt_lines))

    print(f"File saved as: {output_path}")

modify_txt_file(homeDB_path, total_data)

