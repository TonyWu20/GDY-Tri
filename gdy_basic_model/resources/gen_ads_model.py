"""
Prepare molecule insertion information to write MS perl script
"""
from pathlib import Path
from p_tqdm import p_map
from fire import Fire
from graphdiyne import msi_assemble


def main(use_mol: str, lat_dir_pattern: str, mol_z: float = 1.54221):
    """
    Main program to execute
    Args:
        use_mol (str): molecule to be used
        lat_dir_pattern (str): lattice file directory pattern
        mol_z (float): molecule z coord
    """
    lattice_dir = list(Path("msi_models/").rglob(lat_dir_pattern))
    mol_dir = Path("graphdiyne/molecules_models")
    mol_files = list(mol_dir.rglob(use_mol))
    sites = ['metal', 'c1', 'c2', 'c3', 'c4', 'c5']
    execute_arrays = [(mol, lat, site) for mol in mol_files for site in sites
                      for lat in lattice_dir]

    # / Generate all possible combinations among adsorbates, base lattices and
    # / adsorption sites /

    def gen_msi(conf):
        mol_file, lat_dir, site = conf
        #Unpack variable. Single variable for ease of calling p_map()
        factory = msi_assemble.ModelFactory(mol_file, mol_z, lat_dir, site)
        files = list(lat_dir.rglob('*.msi'))
        for file in files:
            factory.assemble_mol(file)

    p_map(gen_msi, execute_arrays)
    print("Done")


def all_molecules(lat_dir_pattern: str, mol_z: float = 1.54221):
    """
    Generate .msi files for all molecules
    """
    main("*.msi", lat_dir_pattern, mol_z)


def to_xsd_scripts(dir_pattern: str):
    """
    Generate msi-to-xsd convert scripts.
    Args:
        dir_pattern (str): directory pattern to match
    """
    root = Path.cwd()
    files = list(root.glob(f"{dir_pattern}/**/*.msi"))
    truncated_paths = [item.relative_to(root) for item in files]
    headlines = ('#!perl\n'
                 'use strict;\n'
                 'use Getopt::Long;\n'
                 'use MaterialsScript qw(:all);\n')
    path_strings = [
        f'"{str(item.parent/item.stem)}"' for item in truncated_paths
    ]
    files_text = ", ".join(path_strings)
    array_text = ("my @params = (\n" + files_text + ');\n')
    actions = ("foreach my $item (@params) {"
               '    my $doc = $Documents{"${item}.msi"};\n'
               "    $doc->CalculateBonds;\n"
               '    $doc->Export("${item}.xsd");\n'
               "    $doc->Save;\n"
               "    $doc->Close;\n"
               "}")
    contents = headlines + array_text + actions
    with open('msi_to_xsd.pl', 'w', newline='\r\n') as file:
        file.write(contents)
    print("Done!")


main_dict = {"model": all_molecules, "to_xsd": to_xsd_scripts}

if __name__ == "__main__":
    Fire(main_dict)
