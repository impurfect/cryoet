"""Download EMPIAR-10491 (5 tilt series) and the EMD-15854 template map."""
import subprocess
from config import DATA, SERIES

FTP = "ftp://ftp.ebi.ac.uk/empiar/world_availability/10491/data"
EMDB = "https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-15854/map/emd_15854.map.gz"


def wget(url, dest, listfile=None):
    dest.mkdir(parents=True, exist_ok=True)
    cmd = ["wget", "-N", "-q", "--show-progress", "-P", str(dest)]
    cmd += ["-i", str(listfile)] if listfile else [url]
    subprocess.run(cmd, check=False)


print("start")

wget(f"{FTP}/gain_ref.mrc", DATA)
print("gain reference")

for s in SERIES:
    wget(f"{FTP}/tiltseries/mdoc/{s}.mrc.mdoc", DATA / "mdoc")
print(f"{len(SERIES)} mdoc files")

# The mdoc files name every movie needed; take the list from them.
names = set()
for m in (DATA / "mdoc").glob("*.mdoc"):
    for line in m.read_text().splitlines():
        if line.strip().startswith("SubFramePath"):
            names.add(line.split("=", 1)[1].strip().replace("\\", "/").rsplit("/", 1)[-1])

urls = DATA / "_urls.txt"
urls.write_text("\n".join(f"{FTP}/tiltseries/data/{n}" for n in sorted(names)) + "\n")
wget(None, DATA / "frames", listfile=urls)
urls.unlink()
print(f"{len(names)} movies")

wget(EMDB, DATA)
subprocess.run(["gunzip", "-f", str(DATA / "emd_15854.map.gz")], check=False)
print("template map")

print(f"done -> {DATA}")
