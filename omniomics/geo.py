"""GEO data access — download processed supplementary files by accession."""
import os, subprocess

def _series_prefix(acc: str) -> str:
    # GSE57575 -> GSE57nnn ; GSE77003 -> GSE77nnn
    num = acc[3:]
    return "GSE" + num[:-3] + "nnn"

def supp_url(acc: str, filename: str | None = None) -> str:
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{_series_prefix(acc)}/{acc}/suppl/"
    return base + filename if filename else base

def list_supp(acc: str) -> list[str]:
    import urllib.request, re
    html = urllib.request.urlopen(supp_url(acc), timeout=60).read().decode("utf-8", "ignore")
    files = re.findall(r'href="([^"]+)"', html)
    return [f for f in files if f.startswith(acc)]

def download(acc: str, filename: str, dest_dir: str) -> str:
    os.makedirs(dest_dir, exist_ok=True)
    out = os.path.join(dest_dir, filename)
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        subprocess.run(["curl", "-sL", supp_url(acc, filename), "-o", out], check=True)
    return out

def download_raw_tar(acc: str, dest_dir: str) -> str:
    return download(acc, f"{acc}_RAW.tar", dest_dir)
