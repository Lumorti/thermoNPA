
# Python script to open a pkl file and print its contents
import pickle
import sys
import os

def open_pkl_file(file_path):
    """
    Open a pickle file and print its contents.
    
    Args:
        file_path (str): Path to the pickle file.
    
    Returns:
        None
    """
    if not os.path.exists(file_path):
        print(f"File {file_path} does not exist.")
        return

    with open(file_path, 'rb') as file:
        data = pickle.load(file)
        print(data)
        print("Data loaded successfully.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python open.py <path_to_pkl_file>")
        sys.exit(1)

    file_path = sys.argv[1]
    open_pkl_file(file_path)
