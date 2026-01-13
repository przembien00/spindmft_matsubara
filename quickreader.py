import h5py
import sys

def read_hdf5_file(file_path,print_dsets=False):
    try:
        with h5py.File(file_path, 'r') as file:
            def print_attrs(name, obj):
                if isinstance(obj, h5py.Dataset):
                    if print_dsets:
                        print(f"  Dataset: {name} -> {obj[...]}")
                    else:
                        print(f"  Dataset: {name} -> shape: {obj.shape}, dtype: {obj.dtype}")
                elif isinstance(obj, h5py.Group):
                    print(f"\033[1mGroup: {name}\033[0m")
                for key, val in obj.attrs.items():
                    print(f"    Attribute: {key:<30} -> {val}")
            file.visititems(print_attrs)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")



if len(sys.argv) > 1:
    file_path = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2] == "showdsets":
        read_hdf5_file(file_path,print_dsets=True)
    else:
        read_hdf5_file(file_path)
else:
    print("to read a hdf file, run: python quickreader.py <filename.hdf5>")
