"""One palette for every figure: satellite classes, DYZ1:DYZ2 ratio and assembly sources (as in Extended Data Fig. 1)."""
ORDER = ['DYZ3', 'DYZ17', 'DYZ19', 'Yq12', 'DYZ18S', 'DYZ1', 'DYZ2', 'DYZ1:DYZ2']   # along chrY; ratio after its parts, in every figure
CLASS = {'DYZ3': '#1f4e79', 'DYZ17': '#8e44ad', 'DYZ19': '#c0392b', 'Yq12': '#333333', 'DYZ18S': '#e377c2', 'DYZ1': '#e6a23c', 'DYZ2': '#2e8b57'}
RATIO = '#8c564b'   # DYZ1:DYZ2
SOURCE = {'Hallast': '#17becf', 'HPRC R2': '#bcbd22'}
KEY = {k.lower(): v for k, v in CLASS.items()}   # dyz3, yq12, ... -> colour
TITLE = lambda size: {'weight': 'bold', 'size': size}   # legend titles: bold, same size as labels
NAMES = {'class': 'Classes', 'source': 'Pangenome', 'lineage': 'Y lineages', 'annotation': 'Annotation'}   # legend titles in every figure
import matplotlib as _mpl
MM = 1 / 25.4   # inches per mm
FS = {'tick': 6.5, 'label': 7, 'legend': 7, 'title': 7, 'point': 6.5, 'letter': 8}   # pt per text role, every figure
_mpl.rcParams.update({'legend.frameon': False, 'legend.handlelength': 1.0, 'legend.handleheight': 0.7, 'legend.handletextpad': 0.4, 'legend.columnspacing': 0.8,   # one legend style
                      'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'font.size': FS['point'], 'pdf.fonttype': 42, 'ps.fonttype': 42,
                      'axes.titlesize': FS['title'], 'axes.labelsize': FS['label'], 'xtick.labelsize': FS['tick'], 'ytick.labelsize': FS['tick'], 'legend.fontsize': FS['legend'], 'legend.title_fontsize': FS['legend'],
                      'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
                      'xtick.major.size': 2, 'ytick.major.size': 2, 'xtick.minor.size': 1.2, 'ytick.minor.size': 1.2, 'xtick.minor.width': 0.4, 'ytick.minor.width': 0.4, 'xtick.major.pad': 2, 'ytick.major.pad': 2, 'axes.labelpad': 2, 'axes.titlepad': 3,
                      'lines.linewidth': 0.75, 'grid.linewidth': 0.4, 'patch.linewidth': 0.5, 'hatch.linewidth': 0.4, 'savefig.dpi': 600, 'savefig.facecolor': 'white'})   # journal style: 5-7 pt, 0.25-1 pt lines


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
    (pdf.savefig(fig) if pdf else fig.savefig(path))
    print(f'written {path} {w / MM:.1f} x {h / MM:.1f} mm' + (f'; WARNING content outside the page (mm): {over}; texts: {out[:4]}' if over else ''), flush=True)
