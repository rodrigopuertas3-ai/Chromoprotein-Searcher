import re
import json
import requests
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

UNIPROT_JSON = "https://rest.uniprot.org/uniprotkb/{acc}.json"
XLSX_NAME = "cromoproteinas.xlsx"

def extract_accession(url: str) -> Optional[str]:
    m = re.search(r"/uniprotkb/([A-Z0-9]+)", url)
    if m:
        return m.group(1)
    if re.fullmatch(r"[A-Z0-9]{6,10}", url.strip()):
        return url.strip()
    return None

def fetch_uniprot_json(accession: str) -> Optional[Dict[str, Any]]:
    url = UNIPROT_JSON.format(acc=accession)
    try:
        r = requests.get(url, timeout=25)
        r.raise_for_status()
        return r.json()
    except Exception:
        return None

NBSP_CHARS = (
    "\u00A0"  
    "\u2007"  
    "\u202F"  
)

def normalize_text(s: str) -> str:
    if not s:
        return ""
    s = re.sub(r"<[^>]+>", " ", s)
    for ch in NBSP_CHARS:
        s = s.replace(ch, " ")
    s = re.sub(r"\s+", " ", s)
    return s.strip()

ABS_STRICT_RE = re.compile(
    r"Abs\s*\(\s*max\s*\)\s*=\s*(\d{2,4})\s*nm",
    flags=re.I
)

EM_STRICT_RE = re.compile(
    r"(?:fluorescence\s+)?emission\s+spectrum\s+which\s+peaks?\s+at\s*(\d{2,4})\s*nm",
    flags=re.I
)

def find_em_strict(entry: Dict[str, Any]) -> str:
    for c in entry.get("comments", []) or []:
        for t in c.get("texts", []) or []:
            val = normalize_text(t.get("value", "") or "")
            m = EM_STRICT_RE.search(val)
            if m:
                return f"{m.group(1)} nm"
    blob = normalize_text(json.dumps(entry, ensure_ascii=False))
    m = EM_STRICT_RE.search(blob)
    if m:
        return f"{m.group(1)} nm"
    return "No encontrado"

def find_abs_in_html(accession: str) -> str:
    import re
    import requests

    NBSP_CHARS = ("\u00A0", "\u2007", "\u202F")
    def _norm(s: str) -> str:
        if not s:
            return ""
        s = re.sub(r"<[^>]+>", " ", s)          
        for ch in NBSP_CHARS:
            s = s.replace(ch, " ")              
        s = re.sub(r"\s+", " ", s).strip()      
        return s

    url = f"https://www.uniprot.org/uniprotkb/{accession}/entry"
    try:
        r = requests.get(url, timeout=25, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        html = r.text
    except Exception:
        return "No encontrado"

    text = _norm(html)

    # intenta anclar al bloque de 'Absorption' y buscar el Abs(max) ahí cerca
    block_pat = re.compile(r"Absorption\s*[:\-]?\s*(.{0,300})", re.I)  
    mblock = block_pat.search(text)
    if mblock:
        window = mblock.group(0)  
        m = re.search(r"Abs\s*\(\s*max\s*\)\s*=\s*(\d{2,4})\s*nm\b", window, re.I)
        if m:
            return f"{m.group(1)} nm"

    # busca el patrón en toda la página (por si el layout varía)
    m = re.search(r"Abs\s*\(\s*max\s*\)\s*=\s*(\d{2,4})\s*nm\b", text, re.I)
    if m:
        return f"{m.group(1)} nm"

    return "No encontrado"

###

def get_ids(entry: Dict[str, Any]) -> Tuple[str, List[str], List[str]]:
    """
    Devuelve:
      - accession (UniProt)
      - embl_like_ids: IDs que sirven para ENA/EMBL/GenBank (DB='EMBL' o 'DDBJ')
      - refseq_like_ids: IDs RefSeq/GenBank (DB='RefSeq' o 'GenBank')
    """
    acc = entry.get("primaryAccession") or ""
    embl_like, refseq_like = [], []
    for xref in entry.get("uniProtKBCrossReferences", []) or []:
        db = (xref.get("database") or "").upper()
        xid = xref.get("id") or ""
        if not xid:
            continue
        if db in {"EMBL", "DDBJ"}:
            if xid not in embl_like:
                embl_like.append(xid)
        elif db in {"REFSEQ", "GENBANK"}:
            if xid not in refseq_like:
                refseq_like.append(xid)
    return acc, embl_like, refseq_like

def get_names(entry: Dict[str, Any]) -> Tuple[str, str]:
    pdsc = entry.get("proteinDescription") or {}
    rec = pdsc.get("recommendedName") or {}
    full = (rec.get("fullName") or {}).get("value", "") if isinstance(rec, dict) else ""
    short = ""
    if isinstance(rec, dict):
        snames = rec.get("shortNames") or []
        if snames and isinstance(snames[0], dict):
            short = snames[0].get("value", "")
    if not short:
        for alt in pdsc.get("alternativeNames") or []:
            snames = alt.get("shortNames") or []
            if snames and isinstance(snames[0], dict):
                short = snames[0].get("value", "")
                break
    return full, short

def get_sequence_info(entry: Dict[str, Any]) -> Tuple[str, int, Optional[int]]:
    seq = entry.get("sequence") or {}
    aa = seq.get("value", "") or ""
    length = seq.get("length", 0) or 0
    mass = seq.get("mass")
    return aa, length, mass

def get_organism(entry: Dict[str, Any]) -> str:
    return (entry.get("organism") or {}).get("scientificName", "") or ""

### PI Teorico

PK = {"Cterm": 2.34,"Nterm": 9.69,"C": 8.33,"D": 3.86,"E": 4.25,"H": 6.00,"K": 10.50,"R": 12.00,"Y": 10.07}

def net_charge(seq: str, pH: float) -> float:
    from collections import Counter
    comp = Counter(seq.upper())
    nterm = 1.0/(1+10**(pH-PK["Nterm"]))
    K = comp["K"]/(1+10**(pH-PK["K"]))
    R = comp["R"]/(1+10**(pH-PK["R"]))
    H = comp["H"]/(1+10**(pH-PK["H"]))
    cterm = -1.0/(1+10**(PK["Cterm"]-pH))
    D = -comp["D"]/(1+10**(PK["D"]-pH))
    E = -comp["E"]/(1+10**(PK["E"]-pH))
    C = -comp["C"]/(1+10**(PK["C"]-pH))
    Y = -comp["Y"]/(1+10**(PK["Y"]-pH))
    return nterm+K+R+H+cterm+D+E+C+Y

def calc_pI(seq: str) -> Optional[float]:
    if not seq: return None
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid=(lo+hi)/2
        charge=net_charge(seq,mid)
        if abs(charge)<1e-4: return round(mid,4)
        if charge>0: lo=mid
        else: hi=mid
    return round((lo+hi)/2,4)

### Referencias

def ref_title_and_link(entry: Dict[str, Any]) -> Tuple[str,str]:
    refs=entry.get("references",[]) or []
    if not refs: return "",""
    cit=refs[0].get("citation",{}) or {}
    title=cit.get("title") or ""
    url=""
    doi=cit.get("doi")
    if doi: url=f"https://doi.org/{doi}"
    else:
        for cr in cit.get("crossReferences",[]) or []:
            if cr.get("database")=="PubMed" and cr.get("id"):
                url=f"https://pubmed.ncbi.nlm.nih.gov/{cr['id']}/"; break
    return title,url

def build_references(entry: Dict[str, Any], max_n:int=3)->str:
    refs=entry.get("references",[]) or []
    out=[]
    for ref in refs[:max_n]:
        cit=ref.get("citation",{}) or {}
        title=cit.get("title") or ""
        journal=cit.get("journal") or ""
        pubdate=cit.get("publicationDate") or ""
        year=re.match(r"(\d{4})",pubdate or "")
        doi=cit.get("doi")
        parts=[title]
        if journal or year: parts.append(f"{journal} ({year.group(1) if year else ''})")
        if doi: parts.append(f"DOI:{doi}")
        s=", ".join([p for p in parts if p])
        if s: out.append(s)
    return " | ".join(out)

### ADN

def _fasta_to_seq(txt: str) -> str:
    lines = [ln.strip() for ln in txt.splitlines() if ln.strip()]
    if not lines: return ""
    if lines[0].startswith(">"): lines = lines[1:]
    return "".join(lines).upper()

def fetch_dna_ncbi(ids: List[str]) -> str:
    """Prueba NCBI (RefSeq/GenBank) con cada ID."""
    for _id in ids:
        try:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            params = {"db":"nuccore","id":_id,"rettype":"fasta","retmode":"text"}
            r = requests.get(url, params=params, timeout=25)
            if r.ok and ">" in r.text:
                seq = _fasta_to_seq(r.text)
                if seq:
                    return seq
        except Exception:
            pass
    return ""

def fetch_dna_ena(ids: List[str]) -> str:
    """Prueba ENA/EMBL con cada ID."""
    for _id in ids:
        try:
            url = f"https://www.ebi.ac.uk/ena/browser/api/fasta/{_id}?download=false"
            r = requests.get(url, timeout=25)
            if r.ok and r.text.startswith(">"):
                seq = _fasta_to_seq(r.text)
                if seq:
                    return seq
        except Exception:
            pass
    return ""

def fetch_best_dna(embl_ids: List[str], refseq_like_ids: List[str]) -> str:
    # 1) NCBI
    if refseq_like_ids:
        seq = fetch_dna_ncbi(refseq_like_ids)
        if seq:
            return seq
    # 2) ENA/EMBL
    if embl_ids:
        seq = fetch_dna_ena(embl_ids)
        if seq:
            return seq
    return ""

### nuevo excel

def save_as_new_sheet(df: pd.DataFrame, filename:str):
    sheet="run_"+datetime.now().strftime("%Y%m%d_%H%M")
    sheet=sheet[:31]
    from openpyxl import load_workbook
    fp=Path(filename)
    if fp.exists():
        with pd.ExcelWriter(filename,engine="openpyxl",mode="a",if_sheet_exists="new") as w:
            df.to_excel(w,index=False,sheet_name=sheet)
    else:
        with pd.ExcelWriter(filename,engine="openpyxl") as w:
            df.to_excel(w,index=False,sheet_name=sheet)
    print(f"✅ Datos guardados en {filename} (hoja: {sheet})")

def fetch_uniprot_txt(accession: str) -> str:
    """
    Descarga el flatfile TXT de UniProt para un accession dado.
    Ej: https://rest.uniprot.org/uniprotkb/P42212.txt
    """
    import requests
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.txt"
    try:
        r = requests.get(url, timeout=20, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        return r.text
    except Exception:
        return ""


def find_abs_from_txt(accession: str) -> str:
    """
    Lee el flatfile TXT de UniProt y devuelve SOLO 'Abs(max) = ### nm'
    desde el bloque 'ABSORPTION'. Si no aparece exactamente así -> 'No encontrado'.
    """
    import re, requests

    url = f"https://rest.uniprot.org/uniprotkb/{accession}.txt"
    try:
        r = requests.get(url, timeout=20, headers={"User-Agent": "Mozilla/5.0"})
        r.raise_for_status()
        txt = r.text
    except Exception:
        return "No encontrado"

    # juntar las lineas
    cc_lines = []
    for ln in txt.splitlines():
        if ln.startswith("CC"):
            cc_lines.append(ln[2:].strip())  # quita 'CC'
    cc_blob = " ".join(cc_lines)
    cc_blob = re.sub(r"\s+", " ", cc_blob)

    mblock = re.search(r"-!-\s*ABSORPTION\s*:\s*(.*?)(?=-!-\s*[A-Z]+:|$)", cc_blob, re.I)
    block = mblock.group(1) if mblock else ""

    m = re.search(r"Abs\s*\(\s*max\s*\)\s*=\s*(\d{2,4})\s*nm\b", block, re.I)
    if m:
        return f"{m.group(1)} nm"

    return "No encontrado"

def find_abs_all(entry: Dict[str, Any]) -> List[str]:
    found = []
    # ABSORPTION estructurado
    for c in entry.get("comments", []) or []:
        if c.get("type") == "ABSORPTION":
            abs_block = c.get("absorption") or {}
            mx = abs_block.get("max")
            if mx:
                s = str(mx).strip()
                if s and s not in found:
                    # Asegura formato "### nm"
                    if not s.lower().endswith("nm"):
                        s = s + " nm"
                    found.append(s)
            # textos dentro
            for t in c.get("texts", []) or []:
                vals = _scan_text_for_abs_all(t.get("value", "") or "")
                for v in vals:
                    if v not in found:
                        found.append(v)
    # otros comentarios
    for c in entry.get("comments", []) or []:
        for t in c.get("texts", []) or []:
            vals = _scan_text_for_abs_all(t.get("value", "") or "")
            for v in vals:
                if v not in found:
                    found.append(v)
    # Global
    if not found:
        blob = json.dumps(entry, ensure_ascii=False)
        for v in _scan_text_for_abs_all(blob):
            if v not in found:
                found.append(v)
    return found

### function main

def main():
    print("🔬 chromoprotein_searcher 🔬")
    links_input=input("👉 Ingresa uno o varios links de UniProt (separados por comas): ").strip()
    links=[l.strip() for l in links_input.split(",") if l.strip()]

    rows=[]
    for link in links:
        acc=extract_accession(link)
        entry=fetch_uniprot_json(acc) if acc else None
        if not entry:
            rows.append({
                "ID de secuencia": acc or "",
                "EMBL/GenBank ID(s)": "",
                "Nombre cromoproteína": "",
                "Nombre corto": "",
                "Secuencia (aa)": "",
                "Longitud secuencia (aa)": "",
                "pI (teórico)": "",
                "Secuencia ADN": "",
                "Longitud secuencia ADN (nt)": "",
                "Organismo de origen": "",
                "Absorbancia (nm)": "No encontrado",
                "Emisión (nm)": "No encontrado",
                "Estudio (título con enlace)": "",
                "Referencias (hasta 3)": "",
                "Peso molecular (Da)": "",
                "Longitud (aa) UniProt": "",
                "Link UniProt": link
            })
            continue

        accession, embl_ids, refseq_ids = get_ids(entry)
        name_full,name_short = get_names(entry)
        aa_seq,length,mass = get_sequence_info(entry)
        org = get_organism(entry)

        abs_val = find_abs_from_txt(accession)

        refs = build_references(entry)
        pI_val = calc_pI(aa_seq)

        # best fetch DNA
        dna_seq = fetch_best_dna(embl_ids, refseq_ids)
        dna_len = len(dna_seq) if dna_seq else ""

        # referencia principal
        title,url = ref_title_and_link(entry)
        if title and url:
            safe_title = title.replace('"', "'")
            study = f'=HYPERLINK("{url}","{safe_title}")'
        else:
            study = title

        rows.append({
            "ID de secuencia": accession,
            "EMBL/GenBank ID(s)": ";".join(embl_ids + refseq_ids) if (embl_ids or refseq_ids) else "",
            "Nombre cromoproteína": name_full,
            "Nombre corto": name_short,
            "Secuencia (aa)": aa_seq,
            "Longitud secuencia (aa)": len(aa_seq),
            "pI (teórico)": pI_val if pI_val is not None else "",
            "Secuencia ADN": dna_seq,
            "Longitud secuencia ADN (nt)": dna_len,
            "Organismo de origen": org,
            "Absorbancia (nm)": abs_val,
            "Estudio (título con enlace)": study,
            "Referencias (hasta 3)": refs,
            "Peso molecular (Da)": mass if mass else "",
            "Longitud (aa) UniProt": length,
            "Link UniProt": f"https://www.uniprot.org/uniprotkb/{accession}/entry"
        })

    df=pd.DataFrame(rows)
    save_as_new_sheet(df,XLSX_NAME)

if __name__=="__main__":
    main()

