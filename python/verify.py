"""Independent checks of the ADS results for Homework 3.

    python verify.py            # prints the comparison table, writes ../docs/py_*.png

Q1  BJT large-signal transconductance and harmonic distortion (Bessel theory)
Q2  emitter-degenerated differential pair: DC transfer curves
Q5  diode envelope detector: ripple vs. RC time constant
"""

import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.integrate import solve_ivp  # noqa: E402
from scipy.optimize import brentq  # noqa: E402
from scipy.special import iv  # noqa: E402

DOCS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs")
VT = 1.380649e-23 * 298.15 / 1.602176634e-19             # ADS default temperature 25 C
rows = []


def report(q, quantity, ads, python):
    rows.append(f"| {q} | {quantity} | {ads} | {python} |")


# --------------------------------------------------------------------- Q1
def q1_bessel():
    """i_C = I_Q exp(x cos wt) = I_Q [I0(x) + 2 sum In(x) cos(n wt)],  x = v_in / V_T."""
    def thd(x, n_max=9):
        return np.sqrt(sum(iv(n, x) ** 2 for n in range(2, n_max + 1))) / iv(1, x)

    x1 = 1e-3 / VT
    report("Q1", "THD for v_in = 1 mV (x = v_in/V_T)", "“97.3 %” (text); spectrum: HD2 ≈ −40 dBc",
           f"x = {x1:.4f} → THD = {100*thd(x1):.3f} %  (≈ x/4)")
    report("Q1", "second-harmonic level HD2 = I2(x)/I1(x)", "≈ −40 dBc (spectrum plot)",
           f"{20*math.log10(iv(2, x1)/iv(1, x1)):.1f} dBc")

    v = np.linspace(1e-4, 0.15, 400)
    x = v / VT
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    axs[0].plot(v * 1e3, 2 * iv(1, x) / (x * iv(0, x)), label="fixed DC current (bias held by feedback)")
    axs[0].plot(v * 1e3, 2 * iv(1, x) / x, label="fixed V_BE bias (as in Q1): DC also grows")
    axs[0].set_ylabel("G_m / g_m")
    axs[0].set_yscale("log")
    axs[0].set_title("Large-signal transconductance (ideal exponential)")
    axs[1].plot(v * 1e3, 100 * thd(x))
    axs[1].set_ylabel("THD (%)")
    axs[1].set_title("Total harmonic distortion of the collector current")
    axs[1].set_ylim(0, 120)
    for ax in axs:
        ax.set_xlabel("input amplitude v_in (mV)")
        ax.grid(alpha=0.3, which="both")
    axs[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(DOCS, "py_q1_bessel.png"), dpi=130)
    plt.close(fig)
    x100 = brentq(lambda xx: thd(xx) - 1.0, 1, 20)
    report("Q1", "input where THD reaches 100 % (ideal device)", "≈ 0.17 V (THD peak)",
           f"{x100*VT*1e3:.0f} mV")


# --------------------------------------------------------------------- Q2
def q2_diff_pair():
    iee, rc = 1e-3, 10e3

    def transfer(vid, re):
        # vid = VT ln(I1/I2) + RE (I1 - I2),  I1 + I2 = IEE.  With u = ln(I1/I2):
        # I1 - I2 = IEE tanh(u/2)  ->  vid = VT u + RE IEE tanh(u/2)  (monotonic in u)
        u = brentq(lambda uu: VT * uu + re * iee * math.tanh(uu / 2) - vid, -500, 500)
        return rc * iee * math.tanh(u / 2)                   # Vo = VC2 - VC1

    vid = np.linspace(-1.5, 1.5, 601)
    fig, ax = plt.subplots(figsize=(7, 4))
    for re in (0.0, 250.0, 500.0):
        vo = np.array([transfer(v, re) for v in vid])
        ax.plot(vid, vo, label=f"R_E = {re:.0f} Ω  (I_EE·R_E = {iee*re/VT:.0f} V_T)")
        # linear range: where the small-signal gain is within 10 % of its maximum
        g = np.gradient(vo, vid)
        lin = vid[g > 0.9 * g.max()]
        report("Q2", f"±linear input range, R_E = {re:.0f} Ω (gain within 10 %)", "plot",
               f"±{lin.max()*1e3:.0f} mV, gain {g.max():.0f} V/V")
    ax.set_xlabel("differential input V1 (V)")
    ax.set_ylabel("V_o = V_C2 − V_C1 (V)")
    ax.set_title("Degeneration trades gain for a wider linear range")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
    ax.set_xlim(-1.5, 1.5)
    fig.tight_layout()
    fig.savefig(os.path.join(DOCS, "py_q2_diffpair.png"), dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------- Q5
def q5_envelope_detector():
    """Diode (Is = 1e-14 A) -> R || C, driven by a 1 V, 1 MHz sine."""
    i_s, n, f = 1e-14, 1.0, 1e6

    def simulate(r, c, t_end=10e-6):
        def rhs(t, v):
            vin = math.sin(2 * math.pi * f * t)
            i_d = i_s * (math.exp(min((vin - v[0]) / (n * VT), 60)) - 1)
            return [(i_d - v[0] / r) / c]
        t = np.linspace(0, t_end, 20001)
        sol = solve_ivp(rhs, (0, t_end), [0.0], t_eval=t, method="LSODA", max_step=2e-9)
        return t, sol.y[0]

    t, v = simulate(500, 100e-12)
    late = t > 5e-6
    report("Q5", "ripple max−min, R = 500 Ω, C = 100 pF (τ = 50 ns)", "0.359 V",
           f"{v[late].max()-v[late].min():.3f} V")
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.8))
    axs[0].plot(t * 1e6, np.sin(2 * np.pi * f * t), color="0.7", lw=0.8, label="input")
    for c, lab in ((100e-12, "τ = 50 ns (assignment)"), (10e-9, "τ = 5 µs"), (100e-9, "τ = 50 µs")):
        tt, vv = simulate(500, c)
        axs[0].plot(tt * 1e6, vv, label=lab)
    axs[0].set_xlabel("time (µs)")
    axs[0].set_ylabel("V_out (V)")
    axs[0].set_title("Envelope detector output, 1 MHz carrier")
    axs[0].legend(fontsize=8)
    taus, rip = np.logspace(-8, -4, 13), []
    for tau in taus:
        tt, vv = simulate(500, tau / 500, t_end=max(10e-6, 8 * tau))
        lt = tt > tt[-1] * 0.7
        rip.append(vv[lt].max() - vv[lt].min())
    axs[1].loglog(taus * 1e6, rip, "o-")
    axs[1].axvline(1, color="k", ls=":", lw=0.8)
    axs[1].text(1.1, 0.4, "carrier period", fontsize=8)
    axs[1].set_xlabel("RC (µs)")
    axs[1].set_ylabel("ripple (V)")
    axs[1].set_title("Ripple only drops once RC ≫ 1/f_carrier")
    for ax in axs:
        ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(DOCS, "py_q5_detector.png"), dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    os.makedirs(DOCS, exist_ok=True)
    q1_bessel()
    q2_diff_pair()
    q5_envelope_detector()
    print("| Q | Quantity | ADS (report) | Python |\n|---|---|---|---|")
    print("\n".join(rows))
