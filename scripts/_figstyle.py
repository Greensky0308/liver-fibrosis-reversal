#!/usr/bin/env python3
"""_figstyle.py — 肝纤维化项目出图样式，对齐老大既定 house style
（样式常量集中定义，供各作图脚本复用）
- 字体 Arial；font.size 7、axes/xtick/ytick/legend 6.2
- Okabe-Ito 多彩分类色 CB
- panel 字母 fontsize 8 @ (-0.16, 1.05)；图例 upper right；去外框
- 导出 PDF + TIF(600 dpi, lzw)
"""
import matplotlib as mpl, matplotlib.pyplot as plt

CB = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9', '#999999']
GRAY = '#999999'

def apply_style():
    mpl.rcParams.update({
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica'], 'font.size': 7,
        'axes.linewidth': 0.6, 'axes.labelsize': 6.2, 'xtick.labelsize': 6.2, 'ytick.labelsize': 6.2,
        'legend.fontsize': 6.2, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': False,
        'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
    })

def despine(ax):
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)

def plab(ax, L):
    ax.text(-0.16, 1.05, L, transform=ax.transAxes, fontsize=8, va='bottom', ha='left')

def setx(ax, size=5.8, rot=45, ha='right'):
    ax.tick_params(axis='x', length=2.5)
    for t in ax.get_xticklabels():
        t.set_fontsize(size); t.set_fontweight('normal')
        if rot: t.set_rotation(rot); t.set_ha(ha)

def save(fig, imgdir, name, qa=True):
    if qa:
        import os; qa = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'out', 'figqa')
        os.makedirs(qa, exist_ok=True)
        fig.savefig(os.path.join(qa, f'{name}.png'), dpi=140, bbox_inches='tight')
    fig.savefig(f'{imgdir}/{name}.pdf', dpi=600)
    fig.savefig(f'{imgdir}/{name}.tif', dpi=600, pil_kwargs={'compression': 'tiff_lzw'})
    plt.close(fig); print('saved', name)
