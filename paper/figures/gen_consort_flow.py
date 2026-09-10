"""
CONSORT flow diagram – v3
Fixes:
  • Vertical arrows end at the visual TOP boundary of the box below them
  • More vertical spacing between boxes
  • Horizontal exclusion arrows end at visual LEFT boundary of exclusion box
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

FIGDIR = '/home/user/VessetVectorTest/paper/figures/'

# ── Palette ──────────────────────────────────────────────────────────────────
C_MF = '#EAF4FB'; C_ME = '#2E86AB'   # main box
C_EF = '#FEF3E2'; C_EE = '#E07B39'   # exclusion box
C_FF = '#E8F5E9'; C_FE = '#2E7D32'   # final / analysed box
C_FL = '#333333'                      # flow arrow

# ── Layout ────────────────────────────────────────────────────────────────────
PAD = 0.10          # FancyBboxPatch pad (same everywhere)
MX  = 3.5           # main-box centre x
BW  = 5.8           # main-box width  → nominal right edge at 3.5+2.9 = 6.4
                    #                   visual right edge at 6.4 + 0.10 = 6.50
EX  = 6.70          # exclusion-box nominal left edge (gap 0.20 from visual r-edge)
EW  = 3.10          # exclusion-box width
EX_VIS = EX - PAD  # visual left boundary where arrowhead lands

fig, ax = plt.subplots(figsize=(10, 14))
ax.set_xlim(0, 10)
ax.set_ylim(0, 15)
ax.axis('off')

# ── Primitive helpers ─────────────────────────────────────────────────────────
def draw_main_box(cy, bh, text, fontsize=10, bold=False, fc=C_MF, ec=C_ME):
    """Draw main-flow box; return (visual_top_y, visual_bot_y)."""
    b = FancyBboxPatch((MX - BW/2, cy - bh/2), BW, bh,
                       boxstyle=f'round,pad={PAD}', fc=fc, ec=ec, lw=1.5, zorder=2)
    ax.add_patch(b)
    ax.text(MX, cy, text, ha='center', va='center', fontsize=fontsize,
            fontweight='bold' if bold else 'normal',
            multialignment='center', zorder=3)
    return cy + bh/2 + PAD, cy - bh/2 - PAD   # visual top, visual bottom

def draw_excl_box(cy, eh, text, fontsize=9):
    """Draw exclusion box with nominal left edge at EX."""
    b = FancyBboxPatch((EX, cy - eh/2), EW, eh,
                       boxstyle=f'round,pad={PAD}', fc=C_EF, ec=C_EE, lw=1.2, zorder=2)
    ax.add_patch(b)
    ax.text(EX + EW/2, cy, text, ha='center', va='center',
            fontsize=fontsize, multialignment='center', zorder=3)

def v_arrow(y_from, y_to):
    """Downward arrow from visual-bottom of upper box to visual-top of lower box."""
    ax.annotate('', xy=(MX, y_to), xytext=(MX, y_from),
                arrowprops=dict(arrowstyle='->', color=C_FL,
                                lw=1.5, shrinkA=0, shrinkB=0), zorder=4)

def h_excl_arrow(branch_y):
    """
    Dashed horizontal line from flow (x=MX) to exclusion-box visual left edge.
    Arrowhead lands at EX_VIS with shrinkB=0.
    """
    ax.plot([MX, EX_VIS], [branch_y, branch_y],
            color=C_EE, lw=1.0, ls='--', zorder=1)
    # Short annotate segment to place solid arrowhead exactly at EX_VIS
    ax.annotate('', xy=(EX_VIS, branch_y), xytext=(EX_VIS - 0.45, branch_y),
                arrowprops=dict(arrowstyle='->', color=C_EE,
                                lw=1.5, shrinkA=0, shrinkB=0), zorder=4)

# ── Computed y-positions (bottom-up so everything stacks cleanly) ─────────────
# Null model boxes (bottom)
NULL_CY  = 2.30;  NULL_BH = 1.30
NULL_VT  = NULL_CY + NULL_BH/2 + PAD   # 3.05
NULL_VB  = NULL_CY - NULL_BH/2 - PAD   # 1.55

FAN_Y = 3.35   # horizontal fan bar y (between pattern-box bottom and null-box tops)

# Pattern detection
PAT_CY = 4.15;  PAT_BH = 0.90
PAT_VT, PAT_VB = PAT_CY + PAT_BH/2 + PAD, PAT_CY - PAT_BH/2 - PAD  # 4.70, 3.60

# Analysed (final)
ANA_CY = 5.75;  ANA_BH = 1.00
ANA_VT, ANA_VB = ANA_CY + ANA_BH/2 + PAD, ANA_CY - ANA_BH/2 - PAD  # 6.35, 5.15

# Exclusion branch 2 — midpoint of gap between Inclusion criteria and Analysed
INCL_CY = 7.95;  INCL_BH = 1.10
INCL_VT, INCL_VB = INCL_CY + INCL_BH/2 + PAD, INCL_CY - INCL_BH/2 - PAD  # 8.60, 7.30

EXCL2_Y = (INCL_VB + ANA_VT) / 2   # (7.30 + 6.35)/2 = 6.825

# Enrolled
ENR_CY = 9.75;  ENR_BH = 0.90
ENR_VT, ENR_VB = ENR_CY + ENR_BH/2 + PAD, ENR_CY - ENR_BH/2 - PAD  # 10.30, 9.20

# Exclusion branch 1 — midpoint of gap between Eligibility and Enrolled
ELIG_CY = 12.05;  ELIG_BH = 0.90
ELIG_VT, ELIG_VB = ELIG_CY + ELIG_BH/2 + PAD, ELIG_CY - ELIG_BH/2 - PAD  # 12.60, 11.50

EXCL1_Y = (ELIG_VB + ENR_VT) / 2   # (11.50 + 10.30)/2 = 10.90

# Source dataset (top)
SRC_CY = 13.60;  SRC_BH = 0.90
SRC_VT, SRC_VB = SRC_CY + SRC_BH/2 + PAD, SRC_CY - SRC_BH/2 - PAD  # 14.15, 13.05

# ── Draw boxes and arrows ─────────────────────────────────────────────────────

# 1. Source dataset
draw_main_box(SRC_CY, SRC_BH,
    'Source dataset: Fehring et al. (2013)\nNFP multi-centre cohort (USA)',
    bold=True)

v_arrow(SRC_VB, ELIG_VT)

# 2. Assessed for eligibility
draw_main_box(ELIG_CY, ELIG_BH,
    'Assessed for eligibility\n(regularly cycling women of reproductive age)')

# Exclusion 1
draw_excl_box(EXCL1_Y, 1.00,
    'Excluded:\n• Menopausal / perimenopausal\n• Postpartum\n• Hormonal use',
    fontsize=8.5)
h_excl_arrow(EXCL1_Y)

v_arrow(ELIG_VB, ENR_VT)

# 3. Enrolled (NFP users in source)
draw_main_box(ENR_CY, ENR_BH,
    'NFP users included in source\n(self-selected, cycle-aware women)')

v_arrow(ENR_VB, INCL_VT)

# 4. Applied inclusion criteria
draw_main_box(INCL_CY, INCL_BH,
    'Applied inclusion criteria:\nCycle length ≥ 18 and ≤ 60 days;\nComplete cycle-length records')

# Exclusion 2
draw_excl_box(EXCL2_Y, 0.90,
    'Excluded:\n• Missing cycle-length data\n• Implausible cycle lengths',
    fontsize=8.5)
h_excl_arrow(EXCL2_Y)

v_arrow(INCL_VB, ANA_VT)

# 5. Analysed (final)
draw_main_box(ANA_CY, ANA_BH,
    'Analysed: N = 118 women\n1,554 cycles (mean 13.2, SD 8.7, range 1–45)',
    bold=True, fc=C_FF, ec=C_FE)

v_arrow(ANA_VB, PAT_VT)

# 6. Pattern detection
draw_main_box(PAT_CY, PAT_BH,
    'Pattern detection applied\n(5 cycle-length-based vesset types)')

# Fan: vertical stem from pattern-box bottom to fan bar
ax.plot([MX, MX], [PAT_VB, FAN_Y], color=C_FL, lw=1.5, zorder=1)

# Fan: horizontal bar
NULL_XS = [1.5, 5.0, 8.5]
ax.plot([NULL_XS[0], NULL_XS[-1]], [FAN_Y, FAN_Y], color=C_FL, lw=1.5, zorder=1)

# Null-model boxes
NULL_BW = 2.60
null_cfg = [
    (NULL_XS[0], '#DBEAFE', '#3B82F6', '$H_G$\nGlobal permutation\n(B = 50,000)'),
    (NULL_XS[1], '#F0FFF0', '#22C55E', '$H_{iid}$\nMultinomial null\n(B = 50,000)'),
    (NULL_XS[2], '#FFF3E0', '#F97316', '$H_W$\nWithin-woman perm.\n(B = 50,000)'),
]
for xc, fc, ec, label in null_cfg:
    ax.annotate('', xy=(xc, NULL_VT), xytext=(xc, FAN_Y),
                arrowprops=dict(arrowstyle='->', color=C_FL,
                                lw=1.2, shrinkA=0, shrinkB=0), zorder=4)
    b = FancyBboxPatch((xc - NULL_BW/2, NULL_CY - NULL_BH/2), NULL_BW, NULL_BH,
                       boxstyle=f'round,pad={PAD}', fc=fc, ec=ec, lw=1.2, zorder=2)
    ax.add_patch(b)
    ax.text(xc, NULL_CY, label, ha='center', va='center',
            fontsize=8.5, multialignment='center', zorder=3)

ax.set_title('Participant Flow and Analysis Design',
             fontweight='bold', fontsize=13, pad=12)

plt.tight_layout()
fig.savefig(FIGDIR + 'fig_consort_flow.pdf', bbox_inches='tight')
fig.savefig(FIGDIR + 'fig_consort_flow.png', bbox_inches='tight', dpi=150)
plt.close(fig)
print("Saved.")
