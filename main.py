
import ROOT
from array import array
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

from paa01_v2 import paaFile

from scipy.optimize import curve_fit


# ============================================================
# CONFIGURACIÓN
# ============================================================

DATA_DIR = "data"

SAMPLE_NS = 8

MIN_WIDTH = 4
MIN_GAP = 25

N_BINS = 27


# ============================================================
# ARCHIVOS
# ============================================================

FILES = sorted(Path(DATA_DIR).glob("*.paa"))

if not FILES:
    raise FileNotFoundError(f"No se encontraron archivos .paa en {DATA_DIR}")


# ============================================================
# DETECCIÓN DE PULSOS
# ============================================================

def contar_pulsos(signal, threshold, min_width=MIN_WIDTH, min_gap=MIN_GAP):

    below = np.array(signal) < threshold

    pulsos = 0
    dentro = False
    ancho = 0
    gap = 0

    for val in below:

        if val:

            if not dentro:
                if pulsos == 0 or gap >= min_gap:
                    dentro = True
                    ancho = 1
                else:
                    ancho += 1
            else:
                ancho += 1

            gap = 0

        else:

            if dentro:
                if ancho >= min_width:
                    pulsos += 1

                dentro = False
                ancho = 0

            gap += 1

    if dentro and ancho >= min_width:
        pulsos += 1

    return pulsos


def encontrar_pulsos(signal, threshold, min_width=MIN_WIDTH, min_gap=MIN_GAP):

    below = np.array(signal) < threshold

    regiones = []

    dentro = False
    inicio = None
    ancho = 0
    gap = 0

    for i, val in enumerate(below):

        if val:

            if not dentro:
                if len(regiones) == 0 or gap >= min_gap:
                    dentro = True
                    inicio = i
                    ancho = 1
                else:
                    ancho += 1
            else:
                ancho += 1

            gap = 0

        else:

            if dentro:
                if ancho >= min_width:
                    regiones.append((inicio, i - 1))

                dentro = False
                inicio = None
                ancho = 0

            gap += 1

    if dentro and ancho >= min_width:
        regiones.append((inicio, len(signal) - 1))

    return regiones


# ============================================================
# ANÁLISIS TEMPORAL
# ============================================================

def tiempo_pulso(signal, region):

    ini, fin = region
    segmento = signal[ini:fin + 1]

    idx_local = np.argmin(segmento)

    return ini + idx_local


def delta_t(signal, threshold):

    regiones = encontrar_pulsos(signal, threshold)

    if len(regiones) != 2:
        return None

    t1 = tiempo_pulso(signal, regiones[0])
    t2 = tiempo_pulso(signal, regiones[1])

    return (t2 - t1) * SAMPLE_NS


# ============================================================
# UTILIDADES
# ============================================================

def info_archivo(file_path):

    df = paaFile(str(file_path))

    return {
        "archivo": file_path.name,
        "threshold": df.paaGetThresholdLevel(),
        "pulse_size": df.paaGetPulseSize(),
        "eventos": df.paaGetPulseCount()
    }


def procesar_archivo(file_path):

    data_file = paaFile(str(file_path))

    threshold = data_file.paaGetThresholdLevel()
    n_events = data_file.paaGetPulseCount()

    delta_ts = []
    doble_pulso = []

    for i in range(n_events):

        pulse = data_file.paaGetPulseRP(i)

        if contar_pulsos(pulse, threshold) == 2:

            doble_pulso.append(i)

            dt = delta_t(pulse, threshold)

            if dt is not None:
                delta_ts.append(dt)

    return delta_ts, doble_pulso, data_file, threshold


# ============================================================
# VISUALIZACIÓN
# ============================================================

def plot_fit_histogram(delta_ts, N0_fit, tau_fit, tau_err, pcov, chi2_red,
                       bins=N_BINS,
                       save_path="results/muon_lifetime_fit.png"):

    import os
    os.makedirs("results", exist_ok=True)
    plt.style.use("default")

    fig, ax = plt.subplots(
    figsize=(12, 7),
    facecolor="white"
)

# ax.set_facecolor("white")

    ax.set_facecolor("white")
    # ax_ratio.set_facecolor("white")

    # fig, ax = plt.subplots(figsize=(12, 7))

    # Histograma sin dibujar
    counts, edges = np.histogram(delta_ts, bins=bins)

    centers = (edges[:-1] + edges[1:]) / 2
    errors = np.sqrt(counts)

    # Dibujar barras
    ax.bar(
        centers,
        counts,
        width=np.diff(edges),
        alpha=0.05,
        edgecolor="black",
        align="center",
        label="Datos experimentales"
    )

    # Barras de error
    ax.errorbar(
        centers,
        counts,
        yerr=errors,
        fmt='o',
        capsize=3,
        markersize=4,
        label="Error estadístico"
    )

    # Curva central
    t_fit = np.linspace(0, max(delta_ts), 1000)
    y_fit = exponencial(t_fit, N0_fit, tau_fit)

    # Valor del fit en centros de bin
    # fit_centers = exponencial(centers, N0_fit, tau_fit)

    # Ratio dato/fit
    # ratio = counts / fit_centers

    # # Error propagado
    # ratio_err = errors / fit_centers

    # ==========================
    # Banda usando covarianza
    # ==========================

    # Derivadas parciales
    dN0 = np.exp(-t_fit / tau_fit)

    dtau = (
        N0_fit
        * np.exp(-t_fit / tau_fit)
        * (t_fit / tau_fit**2)
    )

    # Varianza propagada
    var_y = (
        dN0**2 * pcov[0, 0]
        + dtau**2 * pcov[1, 1]
        + 2 * dN0 * dtau * pcov[0, 1]
    )

    sigma_y = np.sqrt(np.maximum(var_y, 0))

    # # Banda 3σ
    # ax.fill_between(
    #     t_fit,
    #     y_fit - 3*sigma_y,
    #     y_fit + 3*sigma_y,
    #     alpha=0.08,
    #     label=r"$3\sigma$"
    # )

    # Banda 2σ
    ax.fill_between(
        t_fit,
        y_fit - 2*sigma_y,
        y_fit + 2*sigma_y,
        alpha=0.14,
        label=r"$2\sigma$"
    )

    # Banda 1σ
    ax.fill_between(
        t_fit,
        y_fit - sigma_y,
        y_fit + sigma_y,
        alpha=0.25,
        label=r"$1\sigma$"
    )

    # Línea central
    ax.plot(
        t_fit,
        y_fit,
        linewidth=2.5,
        label="Ajuste exponencial"
    )

    # Etiquetas
    ax.set_xlabel("Δt (ns)", fontsize=14)
    ax.set_ylabel("Número de eventos", fontsize=14)
    ax.set_title("Distribución experimental de tiempos de decaimiento del muón", fontsize=16)

    # Grid
    ax.grid(True, alpha=0.18)

    # Caja de texto con resultado

    texto = (
    f"$\\tau = {tau_fit/1000:.3f} \\pm {tau_err/1000:.3f}\\ \\mu s$\n"
    f"$\\chi^2/\\nu = {chi2_red:.2f}$\n"
    f"Eventos: {len(delta_ts)}"
)

    ax.text(
        0.76,
        0.5,
        texto,
        transform=ax.transAxes,
        fontsize=13,
        bbox=dict(
        boxstyle="round,pad=0.4",
        facecolor="white",
        edgecolor="gray",
        alpha=0.95)
    )
    


    ax.legend(fontsize=11, frameon=True)

   # ======================================
# PANEL RATIO SIMPLE
# ======================================

    # ax_ratio.errorbar(
    #     centers,
    #     ratio,
    #     yerr=ratio_err,
    #     fmt='o',
    #     markersize=4,
    #     capsize=3
    # )

    # # Línea de referencia
    # ax_ratio.axhline(
    #     1.0,
    #     linestyle='--',
    #     linewidth=1.5
    # )

    # ax_ratio.set_ylabel("Dato/Fit", fontsize=12)
    # ax_ratio.set_xlabel("Δt (ns)", fontsize=14)

    # ax_ratio.grid(True, alpha=0.2)

    # # límites automáticos razonables
    # ax_ratio.set_ylim(
    #     max(0, np.min(ratio - ratio_err) * 0.9),
    #     np.max(ratio + ratio_err) * 1.1
    # )

    # plt.subplots_adjust(hspace=0.05)

    plt.savefig(save_path, dpi=300, bbox_inches="tight")

    print(f"Imagen guardada en: {save_path}")

    plt.show()

# ============================================================
# FIT
# ============================================================

def exponencial(t, N0, tau):
    return N0 * np.exp(-t / tau)

def ajustar_exponencial(delta_ts):

    hist, edges = np.histogram(delta_ts, bins=N_BINS)

    centers = (edges[:-1] + edges[1:]) / 2

    mask = hist > 0

    popt, pcov = curve_fit(
        exponencial,
        centers[mask],
        hist[mask],
        p0=[max(hist), 2200]
    )

    N0_fit, tau_fit = popt
    tau_err = np.sqrt(np.diag(pcov))[1]

    # Chi cuadrado reducido
    y_fit = exponencial(centers[mask], *popt)

    chi2 = np.sum((hist[mask] - y_fit)**2 / hist[mask])

    ndof = len(hist[mask]) - len(popt)

    chi2_red = chi2 / ndof

    return centers, hist, N0_fit, tau_fit, tau_err, pcov, chi2_red



# ============================================================
# RESUMEN DE ARCHIVOS
# ============================================================

print("\n=== CONFIGURACIÓN DE ARCHIVOS ===\n")

for file in FILES:
    print(info_archivo(file))


# ============================================================
# PROCESAMIENTO GLOBAL
# ============================================================

delta_ts_total = []

print("\n=== PROCESANDO ===\n")

for file in FILES:

    print(f"Procesando: {file.name}")

    delta_ts, doble_pulso, data_file, threshold = procesar_archivo(file)

    print(f"   Dobles pulsos: {len(doble_pulso)}")
    print(f"   Δt válidos:    {len(delta_ts)}")

    delta_ts_total.extend(delta_ts)


print("\n=== RESUMEN FINAL ===")
print(f"Total Δt: {len(delta_ts_total)}")


# ============================================================
# RESULTADO
# ============================================================


centers, hist, N0_fit, tau_fit, tau_err, pcov, chi2_red = ajustar_exponencial(delta_ts_total)

# plot_fit_histogram(
#     delta_ts_total,
#     N0_fit,
#     tau_fit,
#     tau_err
# )
plot_fit_histogram(
    delta_ts_total,
    N0_fit,
    tau_fit,
    tau_err,
    pcov,
    chi2_red
)


print(f"\nVida media medida:")
print(f"τ = {tau_fit:.1f} ± {tau_err:.1f} ns")
print(f"τ = {tau_fit/1000:.3f} ± {tau_err/1000:.3f} μs")
print(f'{N0_fit}, {tau_fit}')

