
# %% IMPORTAR LIBRERÍAS

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


# %% GRÁFICA 1: RELAJACIÓN DE C0(t) - METROPOLIS

# Archivo generado por el programa Fortran
archivo = Path("Metropolis_N256_T0.500_MCS50_R200.csv")

# Leer los resultados
datos = pd.read_csv(archivo)

# Extraer las columnas que necesitamos
t = datos["mcs"].to_numpy()
c0 = datos["c0_mc"].to_numpy()
error_c0 = datos["error_c0"].to_numpy()
c0_teorico = datos["c0_teorico"].to_numpy()

# Crear la figura
fig, ax = plt.subplots(figsize=(8, 5))

# Datos de Monte Carlo con barras de error
ax.errorbar(
    t, c0,
    yerr=error_c0,
    fmt="o",
    markersize=3.5,
    capsize=2,
    elinewidth=0.8,
    linestyle="none",
    color="fuchsia",
    ecolor="fuchsia",
    label="Monte Carlo (Metropolis)"
)

# Predicción teórica
ax.plot(
    t, c0_teorico,
    linewidth=1.8,
    color="darkviolet",
    label="Predicción teórica"
)

# Etiquetas y formato
ax.set_xlabel("Tiempo (MCS)", fontsize=12)
ax.set_ylabel(r"$C_0(t)/C_0(0)$", fontsize=12)
ax.set_title(r"Relajación de $C_0(t)$ - Ising 1-facilitado")

ax.grid(True, alpha=0.25)
ax.legend()
fig.tight_layout()

# Guardar la gráfica
fig.savefig("C0_Metropolis_N256_T0.500_MCS500_R200.png", dpi=300, bbox_inches="tight")

# Mostrar la gráfica
plt.show()






# %% GRÁFICA 2: EVOLUCIÓN DE LA ENERGÍA - METROPOLIS

# Archivo de resultados
archivo = Path("Metropolis_N256_T0.500_MCS50_R200.csv")

# Leer los datos
datos = pd.read_csv(archivo)

# Extraer columnas
t = datos["mcs"].to_numpy()
energia = datos["energia_mc"].to_numpy()
error_energia = datos["error_energia"].to_numpy()
energia_eq = datos["energia_equilibrio"].iloc[0]

# Crear figura
fig, ax = plt.subplots(figsize=(8, 5))

# Energía obtenida con Monte Carlo
ax.errorbar(
    t, energia,
    yerr=error_energia,
    fmt="o",
    markersize=3.5,
    capsize=2,
    elinewidth=0.8,
    linestyle="none",
    color="fuchsia",
    ecolor="fuchsia",
    label="Monte Carlo (Metropolis)"
)

# Energía teórica de equilibrio
ax.axhline(
    y=energia_eq,
    color="darkviolet",
    linestyle="--",
    linewidth=1.8,
    label="Energía teórica de equilibrio"
)

# Formato
ax.set_xlabel("Tiempo (MCS)", fontsize=12)
ax.set_ylim(-0.9, -0.65)
ax.set_ylabel(r"Energía por espín $e(t)$", fontsize=12)
ax.set_title("Evolución de la energía - Ising 1-facilitado")

ax.grid(True, alpha=0.25)
ax.legend()

fig.tight_layout()

# Guardar gráfica
fig.savefig(
    "Energia_Metropolis_N256_T0.500_MCS50_R200.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

# %%
