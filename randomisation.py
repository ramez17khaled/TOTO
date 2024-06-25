import os
import random
import tkinter as tk
from tkinter import filedialog, messagebox, simpledialog
import pandas as pd
import json

def shuffle_table_with_new_index(data):
    """
    Shuffle the table and add a new index while preserving the original sample column.

    Args:
    - data (list of lists): The input data representing a table metabolite in columns.

    Returns:
    - list of tuples: Shuffled table with new index added.
      Each tuple contains an index and a row of shuffled data.
    """
    indexed_data = list(enumerate(data))  # Keep original index with the data
    random.shuffle(indexed_data)  # Shuffle the table
    return indexed_data

def recover_original_order(shuffled_data, original_order):
    """
    Recover the original order using the preserved original index.

    Args:
    - shuffled_data (list of tuples): Shuffled data where each tuple contains an index and a row of shuffled data.
    - original_order (dict): Dictionary mapping original indices (as strings) to their respective rows.

    Returns:
    - list of lists: Recovered table in the original order.
    """
    recovered_table = []
    for shuffled_row in shuffled_data:
        original_index = shuffled_row[0]  # Get the original index from the first element
        original_row = original_order[str(original_index)]  # Get original row from original_order dict
        additional_columns = shuffled_row[1][len(original_row):]  # Get any additional columns that were added
        recovered_table.append(original_row + additional_columns)  # Combine original row with additional columns

    return recovered_table

def gui_interaction():
    """
    Handle GUI interaction, shuffle data, save shuffled data,
    and save original order indices in JSON.

    Uses tkinter for GUI dialogs.

    Returns:
    - None
    """
    # Create a tkinter root window
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    try:
        # Ask user for Excel file path
        file_path = filedialog.askopenfilename(title="Select Excel File", filetypes=[("Excel files", "*.xlsx;*.xls")])

        if not file_path:
            messagebox.showerror("Error", "No file selected. Please select an Excel file.")
            return

        # Ask user for sheet name
        sheet_name = simpledialog.askstring("Sheet Name", "Enter the sheet name:")

        if not sheet_name:
            messagebox.showerror("Error", "No sheet name entered. Please enter a sheet name.")
            return

        # Load data from Excel file
        df = pd.read_excel(file_path, sheet_name=sheet_name)

        # Convert dataframe to list of lists (assuming your data is in this format)
        table = df.values.tolist()

        # Option selection dialog
        choice = messagebox.askquestion("Option", "Do you want to shuffle or recover the table?\n\nSelect 'Yes' for shuffle and 'No' for recover.", icon='question')

        if choice == 'yes':
            # Shuffle the table with new index
            shuffled_table_with_index = shuffle_table_with_new_index(table)

            # Extract the shuffled table without original index
            shuffled_table = [row for _, row in shuffled_table_with_index]

            # Save shuffled table as a new Excel file
            save_path = os.path.dirname(file_path)
            original_filename = os.path.basename(file_path)
            shuffled_filename = f"randomised-{original_filename}"
            shuffled_file_path = os.path.join(save_path, shuffled_filename)

            # Add new order column (1, 2, 3, ...)
            for idx, row in enumerate(shuffled_table):
                row.insert(0, idx + 1)  # Insert new order at the beginning of each row

            # Create new DataFrame for shuffled table
            shuffled_df = pd.DataFrame(shuffled_table, columns=["Order"] + df.columns.tolist())
            shuffled_df.to_excel(shuffled_file_path, index=False)  # Save shuffled DataFrame to Excel

            # Save original order indices to a JSON file
            json_filename = f"original_order-{original_filename}.json"
            json_file_path = os.path.join(save_path, json_filename)

            original_order = {str(i): row for i, row in enumerate(table)}
            with open(json_file_path, 'w') as f:
                json.dump(original_order, f, indent=4)

            messagebox.showinfo("Success", f"Shuffled data saved as '{shuffled_filename}' and original order indices saved as '{json_filename}' in the same directory.")

            # Print the shuffled table with new index and original sample order
            print("\nShuffled Table with New Index (Original Sample Order):")
            for row in shuffled_table:
                print(row)

        elif choice == 'no':
            # Ask user for JSON file path for recovery
            json_path = filedialog.askopenfilename(title="Select JSON File", filetypes=[("JSON files", "*.json")])

            if not json_path:
                messagebox.showerror("Error", "No JSON file selected. Please select a JSON file.")
                return

            # Load original order indices from JSON
            with open(json_path, 'r') as f:
                original_order = json.load(f)

            # Recover the original order using the preserved original index
            shuffled_data = [(i, row) for i, row in enumerate(table)]  # Add original indices to the current data
            recovered_table = recover_original_order(shuffled_data, original_order)

            # Save recovered table as a new Excel file
            save_path = os.path.dirname(file_path)
            original_filename = os.path.basename(file_path)
            recovered_filename = f"recovered-{original_filename}"
            recovered_file_path = os.path.join(save_path, recovered_filename)

            # Add new order column (1, 2, 3, ...)
            for idx, row in enumerate(recovered_table):
                row.insert(0, idx + 1)  # Insert new order at the beginning of each row

            # Create new DataFrame for recovered table
            recovered_df = pd.DataFrame(recovered_table, columns=["Order"] + df.columns.tolist())
            recovered_df.to_excel(recovered_file_path, index=False)  # Save recovered DataFrame to Excel

            messagebox.showinfo("Success", f"Recovered data saved as '{recovered_filename}' in the same directory.")

            # Print the recovered table with new order
            print("\nRecovered Table (Original Order with New Index):")
            for row in recovered_table:
                print(row)

        else:
            messagebox.showinfo("Info", "No action selected. Exiting.")

    except Exception as e:
        messagebox.showerror("Error", f"Error loading or processing data:\n{str(e)}")

    # Destroy the tkinter root window after processing
    root.destroy()

def main():
    """
    Main function to initiate GUI interaction for Excel file selection and sheet name input.
    """
    gui_interaction()

if __name__ == "__main__":
    main()
