"""
Author: [Ramiz KHALED] <ramiz.khaled@inserm.fr>
Infrastructure: [MetaToul/MetaboHub]
Purpose: ABN DB creation.

Overview:
---------
This Python script calculates the probabilities of M+1, M+2, and M+3 isotope peaks for a given molecular formula
based on the natural isotopic abundances of elements. It leverages mathematical formulas and combinatorial 
calculations to determine the contributions of heavy isotopes to molecular mass spectrometry analysis.

Key features:
1. **Input Parsing**: The script includes a function to parse molecular formulas (e.g., "C23H45NO4") into 
   element counts.
2. **Isotopic Abundance Database**: The script uses a dictionary to store the natural isotopic abundances for 
   common elements (C, H, N, O).
3. **M+1 Probability Calculation**: Computes the probability of having a single additional heavy isotope 
   (e.g., one ^13C atom or one ^2H atom).
4. **M+2 Probability Calculation**: Computes the probability of two heavy isotopes appearing simultaneously 
   in various combinations, including contributions from ^18O.
5. **M+3 Probability Calculation**: Computes the probability of three heavy isotopes appearing, accounting 
   for multiple combinations of isotopes.
6. **Flexible Input**: The code can be easily adapted to accept input formulas from a file, a user interface, 
   or a database table.
7. **Output**: The probabilities of M+1, M+2, and M+3 peaks are printed for the given molecular formula.

Applications:
- Useful in analytical chemistry, particularly for mass spectrometry and isotope pattern prediction.
- Can be extended to support additional elements or isotopes as needed.

Dependencies:
- `pandas`: For handling data if extended for tabular input (not directly used in this script but easily adaptable).
- `math`: For combinatorial calculations (e.g., `math.comb`) and mathematical operations.
"""

import pandas as pd
from math import comb
import math

isotopic_abundances = {
    'C': {'p': 1.07/100, 'n': 98.891/100},       # C-13 and C-12
    'H': {'p': 0.0156/100, 'n': 99.9844/100},     # H-2 (D) and H-1
    'N': {'p': 0.365/100, 'n': 99.635/100},       # N-15 and N-14
    'O': {'p1': 0.037/100, 'p2': 0.204/100, 'n': 99.759/100}  # O-17, O-18, and O-16
}

def parse_formula(formula):
    """
    Parses a molecular formula string into a dictionary of elements and their counts.

    Args:
    formula (str): Molecular formula string .

    Returns:
    dict: A dictionary with element symbols as keys and atom counts as values.
          Example: "C6H12O6" -> {'C': 6, 'H': 12, 'O': 6}
    """
    import re
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    return {elem: int(count) if count else 1 for elem, count in elements}

def calculate_m_plus_1_probability(formula):
    """
    Calculates the probability of M+1 isotopes for a given molecular formula.

    Args:
    formula (str): Molecular formula string.

    Returns:
    float: Probability of the M+1 isotope.
    """

    composition = parse_formula(formula)

    x = composition.get('C', 0)
    y = composition.get('H', 0)
    w = composition.get('N', 0)
    z = composition.get('O', 0)

    P_C13 = isotopic_abundances['C']['p']
    P_H2 = isotopic_abundances['H']['p']
    P_N15 = isotopic_abundances['N']['p']
    P_O17 = isotopic_abundances['O']['p1']

    # cas: One isotope BUT NOT O18 
    P_M1 = (
        x * P_C13 +  # Contribution from one C-13
        y * P_H2 +   # Contribution from one H-2
        w * P_N15 +  # Contribution from one N-15
        z * P_O17    # Contribution from one O-17
    )

    return P_M1*100


def calculate_m_plus_2_probability(formula):
    """
    Calculates the probability of M+2 isotopes for a given molecular formula.

    Args:
    formula (str): Molecular formula string.

    Returns:
    float: Probability of the M+2 isotope.
    """

    composition = parse_formula(formula)

    x = composition.get('C', 0)
    y = composition.get('H', 0)
    w = composition.get('N', 0)
    z = composition.get('O', 0)

    P_C13 = isotopic_abundances['C']['p']
    P_H2 = isotopic_abundances['H']['p']
    P_N15 = isotopic_abundances['N']['p']
    P_O17 = isotopic_abundances['O']['p1']
    P_O18 = isotopic_abundances['O']['p2']

    P_M2 = 0

    # Case 1: Two heavy isotopes of the same type
    if x >= 2:
        P_M2 += math.comb(x, 2) * (P_C13**2)  # 2 x C-13
    if y >= 2:
        P_M2 += math.comb(y, 2) * (P_H2**2)   # 2 x H-2
    if w >= 2:
        P_M2 += math.comb(w, 2) * (P_N15**2)  # 2 x N-15
    if z >= 2:
        P_M2 += math.comb(z, 2) * (P_O17**2)  # 2 x O-17

    # Case 2: One heavy isotope of one type and one of another
    if x >= 1 and y >= 1:
        P_M2 += x * P_C13 * y * P_H2
    if x >= 1 and w >= 1:
        P_M2 += x * P_C13 * w * P_N15
    if x >= 1 and z >= 1:
        P_M2 += x * P_C13 * z * P_O17
    if y >= 1 and w >= 1:
        P_M2 += y * P_H2 * w * P_N15
    if y >= 1 and z >= 1:
        P_M2 += y * P_H2 * z * P_O17
    if w >= 1 and z >= 1:
        P_M2 += w * P_N15 * z * P_O17

    # Case 3: One O-18 directly contributes to M+2
    if z >= 1:
        P_M2 += z * P_O18

    return P_M2*100


def calculate_m_plus_3_probability(formula):
    """
    Calculates the probability of M+3 isotopes for a given molecular formula.

    Args:
    formula (str): Molecular formula string.

    Returns:
    float: Probability of the M+3 isotope.
    """
    import math

    composition = parse_formula(formula)

    x = composition.get('C', 0)
    y = composition.get('H', 0)
    w = composition.get('N', 0)
    z = composition.get('O', 0)

    P_C13 = isotopic_abundances['C']['p']
    P_H2 = isotopic_abundances['H']['p']
    P_N15 = isotopic_abundances['N']['p']
    P_O17 = isotopic_abundances['O']['p1']
    P_O18 = isotopic_abundances['O']['p2']

    P_M3 = 0

    # Case 1: Three of the same isotope
    if x >= 3:
        P_M3 += math.comb(x, 3) * P_C13**3
    if y >= 3:
        P_M3 += math.comb(y, 3) * P_H2**3
    if w >= 3:
        P_M3 += math.comb(w, 3) * P_N15**3
    if z >= 3:
        P_M3 += math.comb(z, 3) * P_O17**3

    # Case 2: Two of one isotope and one of another
    if x >= 2 and y >= 1:
        P_M3 += math.comb(x, 2) * P_C13**2 * y * P_H2
    if x >= 2 and w >= 1:
        P_M3 += math.comb(x, 2) * P_C13**2 * w * P_N15
    if x >= 2 and z >= 1:
        P_M3 += math.comb(x, 2) * P_C13**2 * z * P_O17
    if y >= 2 and x >= 1:
        P_M3 += math.comb(y, 2) * P_H2**2 * x * P_C13
    if y >= 2 and w >= 1:
        P_M3 += math.comb(y, 2) * P_H2**2 * w * P_N15
    if y >= 2 and z >= 1:
        P_M3 += math.comb(y, 2) * P_H2**2 * z * P_O17
    if w >= 2 and x >= 1:
        P_M3 += math.comb(w, 2) * P_N15**2 * x * P_C13
    if w >= 2 and y >= 1:
        P_M3 += math.comb(w, 2) * P_N15**2 * y * P_H2
    if w >= 2 and z >= 1:
        P_M3 += math.comb(w, 2) * P_N15**2 * z * P_O17
    if z >= 2 and x >= 1:
        P_M3 += math.comb(z, 2) * P_O17**2 * x * P_C13
    if z >= 2 and y >= 1:
        P_M3 += math.comb(z, 2) * P_O17**2 * y * P_H2
    if z >= 2 and w >= 1:
        P_M3 += math.comb(z, 2) * P_O17**2 * w * P_N15

    # Case 3: One of each of three different isotopes
    P_M3 += x * P_C13 * y * P_H2 * w * P_N15
    P_M3 += x * P_C13 * y * P_H2 * z * P_O17
    P_M3 += x * P_C13 * w * P_N15 * z * P_O17
    P_M3 += y * P_H2 * w * P_N15 * z * P_O17

    # Case 4: One O18 paired with one of another isotope
    if x >= 1:
        P_M3 += x * P_C13 * z * P_O18
    if y >= 1:
        P_M3 += y * P_H2 * z * P_O18
    if w >= 1:
        P_M3 += w * P_N15 * z * P_O18
    if z >= 2:
        P_M3 += math.comb(z, 2) * P_O17 * P_O18

    return P_M3*100


# Execution
formula = "C23H45NO4" 
probability_m1 = calculate_m_plus_1_probability(formula)
print(f"The probability of M+1 for {formula} is: {probability_m1:.4f}")

probability = calculate_m_plus_2_probability(formula)
print(f"The probability of M+2 for {formula} is: {probability:.4f}")


probability_m3 = calculate_m_plus_3_probability(formula)
print(f"The probability of M+3 for {formula} is: {probability_m3:.4f}")
