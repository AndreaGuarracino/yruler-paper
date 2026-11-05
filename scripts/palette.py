"""One palette for every figure: satellite classes, DYZ1:DYZ2 ratio and assembly sources (as in Extended Data Fig. 1)."""
ORDER = ['DYZ3', 'DYZ17', 'DYZ19', 'Yq12', 'DYZ18S', 'DYZ1', 'DYZ2', 'DYZ1:DYZ2']   # along chrY; ratio after its parts, in every figure
CLASS = {'DYZ3': '#1f4e79', 'DYZ17': '#8e44ad', 'DYZ19': '#c0392b', 'Yq12': '#333333', 'DYZ18S': '#e377c2', 'DYZ1': '#e6a23c', 'DYZ2': '#2e8b57'}
RATIO = '#8c564b'   # DYZ1:DYZ2
SOURCE = {'Hallast 2026': '#17becf', 'HPRC R2': '#bcbd22'}
SOURCE_LABEL = {'Hallast 2026': 'Y pangenome (Hallast 2026)', 'HPRC R2': 'HPRC R2 (Lucas 2026)'}   # legend text, both bioRxiv preprints
KEY = {k.lower(): v for k, v in CLASS.items()}   # dyz3, yq12, ... -> colour
TITLE = lambda size: {'weight': 'bold', 'size': size}   # legend titles: bold, same size as labels
REF = {'marker': 'D', 's': 5, 'facecolors': 'none', 'linewidths': 0.8, 'alpha': 0.8}   # male with a reference assembly; edge in panel colour
REF_HANDLE = {'marker': 'D', 'markersize': 5 ** 0.5, 'markerfacecolor': 'none', 'markeredgewidth': 0.8, 'linestyle': 'none'}
BATCH = {'Cross batch': 'o', 'Same batch': '^', 'Reverse cross batch': 's'}   # parent-offspring batch groups
NAMES = {'class': 'Classes', 'source': 'Pangenome', 'lineage': 'Y lineages', 'annotation': 'Annotation', 'superpop': 'Superpopulation'}
SUPERPOP = {'AFR': '#E69F00', 'AMR': '#CC79A7', 'EAS': '#56B4E9', 'EUR': '#009E73', 'SAS': '#F0E442'}   # Okabe-Ito, colourblind-safe
SUPERPOP_NAME = {'AFR': 'Africa (AFR)', 'AMR': 'Americas (admixed) (AMR)', 'EAS': 'East Asia (EAS)', 'EUR': 'Europe (EUR)', 'SAS': 'South Asia (SAS)'}   # legend titles in every figure
import os
import matplotlib as _mpl
MM = 1 / 25.4   # inches per mm
PAD = 1.5 * MM   # figure margin, every figure
GRID = '#e8edef'   # grid lines
FS = {'tick': 6, 'label': 7.5, 'legend': 6, 'title': 9, 'point': 6, 'letter': 15}   # Fig. 1 sizes (Alejandro): titles 9 bold, axis labels 7.5 italic, letters 15 bold   # pt per text role, every figure
_mpl.rcParams.update({'legend.frameon': False, 'legend.handlelength': 1.0, 'legend.handleheight': 0.7, 'legend.handletextpad': 0.4, 'legend.columnspacing': 0.8,   # one legend style
                      'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': FS['point'], 'pdf.fonttype': 42, 'ps.fonttype': 42,
                      'axes.titlesize': FS['title'], 'axes.titleweight': 'bold', 'axes.labelsize': FS['label'], 'xtick.labelsize': FS['tick'], 'ytick.labelsize': FS['tick'], 'legend.fontsize': FS['legend'], 'legend.title_fontsize': FS['legend'],
                      'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                      'xtick.major.size': 2, 'ytick.major.size': 2, 'xtick.minor.size': 1.2, 'ytick.minor.size': 1.2, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4, 'xtick.major.pad': 2, 'ytick.major.pad': 2, 'axes.labelpad': 2, 'axes.titlepad': 3,
                      'lines.linewidth': 0.75, 'grid.linewidth': 0.4, 'grid.color': GRID, 'patch.linewidth': 0.5, 'hatch.linewidth': 0.4, 'savefig.dpi': 600, 'savefig.facecolor': 'white'})   # journal style: 5-7 pt, 0.25-1 pt lines


DOSE = '#526774'   # Y dosage
INK, MUTED, LINE, GREY, ORANGE = '#22313a', '#55666f', '#60717a', '#a9b4ba', '#d9622b'   # text, guides, tree lines, other males, assemblies
J2 = '#bd6078'   # J2a1a1 in ED8
BAND = {'DYZ3': '#e8f1f7', 'DYZ19': '#f8e9ec'}   # highlighted rows
LAND = {'face': '#ebeff1', 'edge': '#bbc5ca'}
CONTIG = {'face': '#f0f0f0', 'edge': '#bbbbbb'}   # ED1 contigs
UNITS = {'Other units': '#a584b9', 'Unannotated': '#dbe1e5'}   # ED2 Yq12 composition
SCHEME = {'reference': '#eef3f5', 'train': '#e7f1f7', 'test': '#f9eee7', 'output': '#eaf4ef'}   # Fig. 1a boxes
LINEAGE = {'A': '#a73945', 'B': '#d3b183', 'C': '#2684b2', 'D': '#789bac', 'E': '#de915b', 'F': '#94a6ad', 'G': '#b6a14b', 'H': '#b3a1c8',
           'I': '#8b6baa', 'J': '#439a88', 'L': '#83b955', 'M': '#68a5a7', 'NO': '#b0b6bf', 'N': '#9bbc59', 'O': '#62b9b5', 'Q': '#d8bb51',
           'R': '#77628f', 'S': '#b78d69', 'T': '#c4d39e', 'Unassigned': '#dddddd'}   # Y lineages
HEAT = ['#123c73', '#2a78d6', '#9cc3ee', '#d4d3ce', '#f2a6a5', '#e34948', '#8f1f1f']   # ED3 copy-number classes
STAIN = {'gneg': '#fcfcfb', 'gpos50': '#a9a8a4', 'acen': '#5d5c58', 'gvar': '#d6d5d1'}   # cytoband stains
COLOR = {**CLASS, 'DYZ1:DYZ2': RATIO, 'Y dosage': DOSE}   # every panel title colour


TITLE_PAD = 6.5   # pt between panel title baseline and axes, every figure


def head(ax, text, pad=TITLE_PAD):
    """Panel title, centred, in its class colour."""
    ax.set_title(text, color=COLOR.get(text, 'black'), pad=pad)


def top_legend(fig, x=.5, loc='lower center', gap=1.2, **kw):
    """Figure legend 1.2 mm above the top letter or title."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    T = [t for t in list(fig.texts) + [a.title for a in fig.axes if a.get_visible()] if t.get_visible() and t.get_text().strip()]
    top = max(t.get_window_extent(r).y1 for t in T) / fig.bbox.height
    return fig.legend(loc=loc, bbox_to_anchor=(x, top + gap * MM / fig.get_figheight()), borderaxespad=0, **kw)


def letter(fig, s, x, y):
    """Panel letter at figure fraction x (left), y (top)."""
    fig.text(x, y, s, fontsize=FS['letter'], fontweight='bold', ha='left', va='top')


def save(fig, path, pdf=None):
    """Save at the set page size; report overflow."""
    from matplotlib.transforms import Bbox
    fig.canvas.draw(); r = fig.canvas.get_renderer(); w, h = fig.get_size_inches()
    A = [a for a in fig.axes if a.get_visible()]; A = list(fig.texts) + [t for a in A for t in a.texts + [a.title, a._left_title, a._right_title] + ([a.xaxis.label, a.yaxis.label] if a.axison else [])]
    T = [(t.get_text(), Bbox(t.get_window_extent(r).get_points() / fig.dpi)) for t in A if t.get_visible() and t.get_text().strip()]
    b = Bbox.union([fig.get_tightbbox(r)] + [x for _, x in T])   # tight bbox squeezes axis labels
    out = [s for s, x in T if x.x0 < 0 or x.y0 < 0 or x.x1 > w or x.y1 > h]
    over = {k: round(float(v) / MM, 1) for k, v in (('left', -b.x0), ('bottom', -b.y0), ('right', b.x1 - w), ('top', b.y1 - h)) if v > 0.05 * MM}
    for ax in fig.axes:   # Fig. 1 style: italic axis labels, bold titles
        ax.xaxis.label.set_fontstyle('italic'); ax.yaxis.label.set_fontstyle('italic'); ax.title.set_fontweight('bold')
    (pdf.savefig(fig) if pdf else fig.savefig(path, bbox_inches='tight', pad_inches=PAD))   # same 1.5 mm margin on every side
    if not pdf and str(path).endswith('.pdf'): crop_ink(path)   # margins measured on the ink
    if not pdf and '/main/' in str(path):   # main figures: also native SVG, editable text
        with _mpl.rc_context({'svg.fonttype': 'none'}): fig.savefig(str(path)[:-4] + '.svg'); print(f'written {str(path)[:-4]}.svg', flush=True)
    print(f'written {path} {w / MM:.1f} x {h / MM:.1f} mm' + (f'; WARNING content outside the page (mm): {over}; texts: {out[:4]}' if over else ''), flush=True)


def crop_ink(path, pad=1.5, dpi=300):
    """Crop a PDF to its ink plus pad mm."""
    import subprocess, tempfile, numpy as np
    from PIL import Image
    from pypdf import PdfReader, PdfWriter
    with tempfile.TemporaryDirectory() as t:
        subprocess.run(['pdftoppm', '-r', str(dpi), '-png', '-singlefile', str(path), f'{t}/p'], check=True); a = np.asarray(Image.open(f'{t}/p.png').convert('L'))
    ys, xs = np.nonzero(a < 250); k = 72 / dpi; q = pad / 25.4 * 72
    if not len(xs): raise SystemExit(f'ERROR no ink in {path}')
    r = PdfReader(str(path)); w = PdfWriter(); pg = r.pages[0]; x0, y0, x1, y1 = (float(v) for v in pg.mediabox)
    box = (x0 + xs.min() * k - q, y1 - (ys.max() + 1) * k - q, x0 + (xs.max() + 1) * k + q, y1 - ys.min() * k + q)
    pg.mediabox.lower_left, pg.mediabox.upper_right = box[:2], box[2:]; pg.cropbox.lower_left, pg.cropbox.upper_right = box[:2], box[2:]
    w.add_page(pg); w.write(str(path))
    print(f'cropped {path} to ink + {pad} mm: {(box[2] - box[0]) / 72 * 25.4:.1f} x {(box[3] - box[1]) / 72 * 25.4:.1f} mm', flush=True)
