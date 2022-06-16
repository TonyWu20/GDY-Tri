"""
Write castep running script (lsf) for seed files
Usage:
    python3 write_lsf_scripts.py keyword
"""
from pathlib import Path
from p_tqdm import p_map
from fire import Fire


def get_main_name(dirname: str) -> str:
    """
    Fetch main name, remove suffix like 'opt'
    Arg:
        dirname
    return:
        main_name
    """
    elements = dirname.split("_")
    main_name = str.join("_", elements[:-1])
    return main_name


def write_scripts(dirpath: Path):
    """
    Write lsf running script in the folder
    Arg:
        dir_path (Path): Path of the folder
    """
    names = get_main_name(dirpath.name)
    cmd_prefix = (
        "/home-yw/Soft/msi/MS70/MaterialsStudio7.0/etc/CASTEP/bin/RunCASTEP.sh -np $NP"
    )
    cmd = f"{cmd_prefix} {names}"
    pre = [
        "APP_NAME=intelY_mid",
        "NP=12",
        "NP_PER_NODE=12",
        "OMP_NUM_THREADS=1",
        'RUN="RAW"\n',
    ]
    pre = [f"{item}\n" for item in pre]
    file = f"{str(dirpath)}/MS70_YW_CASTEP.lsf"
    script = open(file, "w")
    for line in pre:
        script.write(line)
    script.write(cmd)
    script.close()


def main(keyword: str):
    """
    Main program
    Arg:
        keyword: Match pattern
    """
    src = Path.cwd()
    dirs = list(src.rglob(keyword))
    print("Writing lsf scripts...")
    p_map(write_scripts, dirs)


if __name__ == "__main__":
    Fire(main)
