# Changelog — padic-balls-animation.html

---

## 2026-10-18 (suite v1.2.0)

- Niveles `n = 0, …, N−1` en `computeAllLevels`, `computeMuFast`, `computeMuFastBit`, el invariante (`μ_E+μ_A+μ_I = N·p^N`), la tabla de perfil, el slider, la tira de historia por nivel, el selector del grafo de bolas y el μ_E del panel *How to*.
- Nivel `n = 0`: el fractal pinta la raíz ℤₚ con su tipo.
- Corrección del texto del centro de una bola: quedan `N−n` genes libres.

---

## 2026-04-27 (sesión Cowork)

### Nuevas funcionalidades

#### 1. Grafo de transición (`renderTransitionGraph`)
- Pestaña **"Transition graph"** en el panel de dinámicas.
- Muestra el grafo funcional $x \to f(x)$ para todos los estados.
- Nodos coloreados por pertenencia a ciclo (cada ciclo tiene su propio color; nodos transitorios en gris).
- Layout: nodos de ciclo en anillo interior, nodos transitorios en corona exterior agrupados por ancestro.
- Etiquetas decimales para redes con $\leq 64$ estados; omitidas para redes más grandes.
- Límite: $\leq 256$ estados (advertencia para redes más grandes).
- Llama automáticamente desde `runCompute()` y desde `switchDynTab('transgraph', ...)`.

#### 2. Estructura polinomial ANF (`renderPolynomialPanel`)
- Pestaña **"Polynomial structure"** en el panel de dinámicas.
- Calcula la **Algebraic Normal Form (ANF)** de cada función de actualización $T_i: \mathbb{F}_2^N \to \mathbb{F}_2$ mediante la transformada de Möbius — $O(N \cdot 2^N)$, en tiempo real para $N \leq 20$.
- **Tabla** con: nombre del gen, polinomio ANF (truncado a 8–12 términos para redes grandes), número de monomios, grado máximo, influence score $I_{\text{poly}}$.
- **Detección automática** de estructuras algebraicas especiales:
  - Identidades: $T_i = x_j$ → genes sincronizados; reducción dimensional en atractores.
  - OR booleano: $a \oplus b \oplus ab = a \vee b$.
  - Factorización con factor lineal común.
- **Nota metodológica** al pie: explica $I_{\text{poly}}$ y la relación con el ordering óptimo.
- **Detección e informe de identidades** ($T_i = x_j, i \neq j$) con comentario sobre reducción dimensional.
- **Grafo de dependencias $G_{\text{poly}}$**: nodo por gen, arista $x_j \to x_i$ si $x_j$ aparece en $T_i$; nodos del gen más influyente en teal, resto en azul/gris.
- Requiere p=2. Para p>2 muestra mensaje informativo.
- Llama automáticamente desde `runCompute()` y desde `switchDynTab('polystruct', ...)`.

#### 3. Soporte para redes custom (CSV) — mejoras previas ya incorporadas
- Selector "Red personalizada (CSV)…" con zona de drag-drop.
- Auto-detección de formato: 2 columnas (enteros decimales) o 2N columnas (vectores de bits little-endian).
- `net.genes` null-safe en tooltip, `renderCyclesPanel`, `renderHowtoPanel`.

---

## Resultados del análisis ANF computacional (2026-04-27)

### Drosophila melanogaster (N=6)
- bcd y tll son identidades puras ($T_0=\text{bcd}$, $T_5=\text{tll}$): morfógenos maternos fijos.
- Simetría bilateral exacta: el intercambio kni↔gt es un automorfismo de $f$.
- Todos los polinomios factorizan con factor lineal (bcd o tll): red descompuesta en dos subsistemas.

### C. elegans (N=8)
- Identidades: $T_{\text{cul-1}} = \text{cdk-2}$ y $T_{\text{fzr-1}} = \text{cdc-14}$. Dimensión efectiva en atractores: 5 (de 8).
- OR booleano: $T_{\text{cki-1}} = \text{fzr-1} \vee \text{cdc-14}$ → colapsa a cdc-14 en atractores.
- cdk-1 tiene 40 monomios de grado $\leq 4$ — integrador central del ciclo.

### Notas para Paper II
- §09d usa genes (Bcd, **Caudal**, Hb, Kr, Kni, Gt) pero el software tiene (bcd, **tll**, hb, Kr, kni, gt). Modelos distintos — hay que corregir §09d.
- §08d incompleto en el paper; los polinomios reales de C. elegans están ahora calculados.
- Ver `notes/anf_analysis_grn_2026-04-27.md` y `notes/ANF-GRN-future-article.md` en la carpeta del coloquio.

---

## Archivos modificados

- `padic-balls-animation.html` — applet principal (transición + ANF + CSV).
- `boolean_N4_example.csv` — red booleana ejemplo N=4 con bolas expansivas.
- `z2mod16_N4_p2.csv` — red z²mod16 en formato vectores de bits.

