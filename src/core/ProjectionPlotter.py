import os
import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math as m

class PlotSpec:
    FIGSIZE_1D = (6, 6)
    FIGSIZE_2D = (7, 7)
    DPI = 300
    MARKER = "o"
    MARKER_SIZE = 20
    COLOR = "tab:red"
    LINEWIDTH = 1.5
#   SCATTER_SIZE_2D = 1
#   SCATTER_EDGE_COLOR_2D = "none"
    CONTOUR_LEVELS = 10
    CMAP = "viridis"
    COLORBAR_SHRINK = 0.82
    GRID_N = 100
    TICK_MARGIN_FRAC = 0.05
    BBOX_INCHES = "tight"
    GRID_METHOD = "linear"
    COLORBAR_ORIENT = "horizontal"
    NT_LONG = 6
    AXIS_LABEL_SIZE = 14
    XLABPAD_1D = 5
    YLABPAD_1D = 10
    XLABPAD_2D = 3
    YLABPAD_2D = 8
    ZLABPAD_2D = 0.12

    GREEK = {
        "alpha": r"\alpha", "beta": r"\beta", "gamma": r"\gamma",
        "delta": r"\delta", "epsilon": r"\epsilon", "zeta": r"\zeta",
        "eta": r"\eta", "theta": r"\theta", "iota": r"\iota",
        "kappa": r"\kappa", "lambda": r"\lambda", "mu": r"\mu",
        "nu": r"\nu", "xi": r"\xi", "omicron": r"o",
        "pi": r"\pi", "rho": r"\rho", "sigma": r"\sigma",
        "tau": r"\tau", "upsilon": r"\upsilon", "phi": r"\phi",
        "chi": r"\chi", "psi": r"\psi", "omega": r"\omega"}

    COST_UNITS = {
        "1": r"$\mathrm{cm}^{-1}$",
        "2": r"$\mathrm{eV}$",
        "3": r"$\mathrm{E_h}$",
        "4": r"$\mathrm{kcal}\,\mathrm{mol}^{-1}$",
        "5": r"$\mathrm{kJ}\,\mathrm{mol}^{-1}$"}

    def __init__(self, inp_dir, filename):
        self.inp_dir = inp_dir
        self.filename = filename
        self._parse_filename()
        self._build_latex()

    def _parse_filename(self):
        name = self.filename[:-4]
        parts = name.split("_")
        self.lab_wcoors = []
        self.unit_wcoors = []
        i = 1
        while i < len(parts) - 2:
            self.lab_wcoors.append(parts[i])
            self.unit_wcoors.append(parts[i + 1])
            i += 2
        self.cost_label = parts[-2]
        self.cost_unit = parts[-1]
        self.ndim = len(self.lab_wcoors)

    def _latexify(self, kind, value):
        if kind == "wcoor_label":
            m0 = re.fullmatch(r"([A-Za-z]+)(\d*)", value)
            if not m0:
                return value
            base, idx = m0.groups()
            base_low = base.lower()
            if base_low in self.GREEK:
                sym = self.GREEK[base_low]
            elif len(base) == 1:
                sym = base
            else:
                sym = f"\\mathrm{{{base}}}"
            if idx:
                return f"${sym}_{{{idx}}}$"
            return f"${sym}$"
        if kind == "wcoor_unit":
            if value == "A":
                return r"$\mathrm{\AA}$"
            if value.lower() == "deg":
                return r"$^\mathrm{o}$"
            return value
        if kind == "cost_label":
            if value.isupper():
                return value
            if re.search(r"[A-Z]", value):
                return " ".join(
                        re.findall(r"[A-Z][^A-Z]*", value))
            return value
        if kind == "cost_unit":
            if self.cost_label.lower() == "energy":
                return self.COST_UNITS.get(value, value)
            return value
        return value

    def _build_latex(self):
        self.latex_wcoors = []
        for l, u in zip(self.lab_wcoors, self.unit_wcoors):
            self.latex_wcoors.append({
                "label": self._latexify("wcoor_label", l),
                "unit": self._latexify("wcoor_unit", u)})
        self.latex_cost = {
            "label": self._latexify("cost_label", self.cost_label),
            "unit": self._latexify("cost_unit", self.cost_unit)}

    def load_data(self):
        path = os.path.join(self.inp_dir, self.filename)
        cost = []
        wcoors = []
        with open(path, "r") as f:
            for line in f:
                line = line.split("|")[0]
                parts = line.split()
                cost.append(float(parts[0]))
                wcoors.append([float(x) for x in parts[1:]])
        self.cost = np.array(cost)
        self.wcoors = np.array(wcoors)

    def _find_ticks(self, Cmin, Cmax, ntarget, fac):
        def ceil_scale(x, S):
            return m.ceil(S * x) / S
        def floor_scale(x, S):
            return m.floor(S * x) / S
        def to_int(S, x):
            return int(S * x)
        def decimals(x):
            s = format(x, "f")
            if "." not in s:
                return 0
            return len(s.rstrip("0").split(".")[1])
        width = Cmax - Cmin
        marg = fac * width
        tol = 2.0 * marg
        cmin0 = Cmin + marg
        cmax0 = Cmax - marg
        d = 0
        sign = 1
        S = 1.0
        a = ceil_scale(cmin0, S)
        b = floor_scale(cmax0, S)
        while True:
            if (sign == -1) ^ (a < b):
                if d == 0:
                    sign = -1
                else:
                    if sign == 1:
                        cmin = a
                        cmax = b
                    break
            cmin = a
            cmax = b
            S *= 10.0 ** sign
            a = ceil_scale(cmin0, S)
            b = floor_scale(cmax0, S)
            if sign == -1 and \
                (a > cmin0 + tol or 
                    b < cmax - tol):
                break
            d += sign
        span = cmax - cmin
        max_dec = max(decimals(cmin), decimals(cmax))
        best_nt = None
        best_ticks = None
        for nt in range(3, ntarget + 1):
            disp = span / (nt - 1)
            ticks = [cmin + i * disp for i in range(nt)]
            ok = True
            for t in ticks:
                if decimals(t) > max_dec + int(nt == 3):
                    ok = False
                    break
            if ok:
                best_nt = nt
                best_ticks = ticks
        best_disp = span / (best_nt - 1)
        best_max_dec = max(max_dec, decimals(best_disp))
        fmt = "{:." + str(best_max_dec) + "f}"
        return [float(fmt.format(t)) for t in best_ticks]

    def plot(self, out_dir):
        if self.ndim == 1:
            idx = np.argsort(self.wcoors[:, 0])
            self.wcoors = self.wcoors[idx]
            self.cost = self.cost[idx]
            self._plot_1d(out_dir)
        elif self.ndim == 2:
            self._plot_2d(out_dir)

    def _plot_1d(self, out_dir):
        x = self.wcoors[:, 0]
        y = self.cost
        xlabel = f"{self.latex_wcoors[0]['label']} / " \
                 f"{self.latex_wcoors[0]['unit']}"
        ylabel = f"{self.latex_cost['label']} / " \
                 f"{self.latex_cost['unit']}"
        plt.figure(figsize=self.FIGSIZE_1D)
        plt.plot(x, y, color=self.COLOR, linewidth=self.LINEWIDTH)
        plt.scatter(x, y, color=self.COLOR, marker=self.MARKER,
                    s=self.MARKER_SIZE)
        ax = plt.gca()
        ax.set_xlabel(xlabel, fontsize=self.AXIS_LABEL_SIZE, 
                      labelpad=self.XLABPAD_1D)
        ax.set_ylabel(ylabel, fontsize=self.AXIS_LABEL_SIZE, 
                      labelpad=self.YLABPAD_1D)
        ax.set_xticks(self._find_ticks(x.min(), x.max(),
                      self.NT_LONG, self.TICK_MARGIN_FRAC))
        ax.set_yticks(self._find_ticks(y.min(), y.max(),
                      self.NT_LONG, self.TICK_MARGIN_FRAC))
        base = os.path.splitext(self.filename)[0]
        out = os.path.join(out_dir, base + ".png")
        plt.savefig(out, dpi=self.DPI, bbox_inches=self.BBOX_INCHES)
        plt.close()

    def _plot_2d(self, out_dir):
        x = self.wcoors[:, 0]
        y = self.wcoors[:, 1]
        z = self.cost
        lx = self.latex_wcoors[0]
        ly = self.latex_wcoors[1]
        ux = self.unit_wcoors[0]
        uy = self.unit_wcoors[1]
        fig, ax = plt.subplots(figsize=self.FIGSIZE_2D)
        x_range = x.max() - x.min()
        y_range = y.max() - y.min()
        if ux == uy:
            if y_range > x_range:
                x, y = y, x
                lx, ly = ly, lx
                ux, uy = uy, ux
                x_range, y_range = y_range, x_range
            ax.set_aspect(y_range / x_range)
            nt_x = self.NT_LONG
            nt_y = int(np.ceil(y_range * self.NT_LONG / x_range))
        else:
            nt_x = self.NT_LONG
            nt_y = self.NT_LONG
        xi = np.linspace(x.min(), x.max(), self.GRID_N)
        yi = np.linspace(y.min(), y.max(), self.GRID_N)
        X, Y = np.meshgrid(xi, yi)
        Z = griddata((x, y), z, (X, Y), method=self.GRID_METHOD)
        Zmax = np.nanmax(z)
        Zfill = Zmax + 0.02 * (Zmax - np.nanmin(z))
        Z = np.where(np.isnan(Z), Zfill, Z)
        cs = ax.contourf(X, Y, Z, levels=self.CONTOUR_LEVELS, 
                        cmap=self.CMAP)
#       ax.scatter(X, Y, c=Z, s=self.SCATTER_SIZE_2D,
#                  edgecolors=self.SCATTER_EDGE_COLOR_2D)
        ax.set_xticks(self._find_ticks(x.min(), x.max(), nt_x,
                      self.TICK_MARGIN_FRAC))
        ax.set_yticks(self._find_ticks(y.min(), y.max(), nt_y,
                      self.TICK_MARGIN_FRAC))
        cbar = plt.colorbar(cs, orientation=self.COLORBAR_ORIENT,
                            pad=self.ZLABPAD_2D,
                            shrink=self.COLORBAR_SHRINK)
        cbar.set_label(f"{self.latex_cost['label']} / "
                       f"{self.latex_cost['unit']}",
                       fontsize=self.AXIS_LABEL_SIZE)
        cbar.ax.set_xticks(self._find_ticks(z.min(), z.max(),
                           self.NT_LONG, self.TICK_MARGIN_FRAC))
        ax.set_xlabel(f"{lx['label']} / {lx['unit']}", 
                      fontsize=self.AXIS_LABEL_SIZE,
                      labelpad=self.XLABPAD_2D)
        ax.set_ylabel(f"{ly['label']} / {ly['unit']}",
                      fontsize=self.AXIS_LABEL_SIZE, 
                      labelpad=self.YLABPAD_2D)
        base = os.path.splitext(self.filename)[0]
        out = os.path.join(out_dir, base + ".png")
        plt.savefig(out, dpi=self.DPI, bbox_inches=self.BBOX_INCHES)
        plt.close()
        print("{:>{w}}".format(out, w=len(out) + 9))

def plotter(inp_dir, out_dir):
    for filename in os.listdir(inp_dir):
        spec = PlotSpec(inp_dir, filename)
        spec.load_data()
        spec.plot(out_dir)

if __name__ == "__main__":
    inp_dir = sys.argv[1]
    out_dir = sys.argv[2]
    plotter(inp_dir, out_dir)
