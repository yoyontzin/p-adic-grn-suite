# Manual de usuario — padic-grn-search

## Instalación

```bash
cd software-suite/padic-grn-search
pip install -e .
```

## Comandos principales

### 1) GA

```bash
pgsearch-ga ../padic-grn-core/examples/transitions_p2_n2_identity.csv 2 2 --objective minimize_mu_E --generations 5 --population 10 --seed 7
```

### 2) B&B mínimo de mu_E

```bash
pgsearch-bb-min-muE ../padic-grn-core/examples/transitions_p2_n2_identity.csv 2 2
```

### 3) B&B máximo de mu_A

```bash
pgsearch-bb-max-muA ../padic-grn-core/examples/transitions_p2_n2_identity.csv 2 2
```

## Recomendaciones

- Para pruebas rápidas usar redes pequeñas (`N<=8`).
- Para certificación exacta usar `pgsearch-bb-*`.
- Para exploración rápida en `N` grande usar `pgsearch-ga` y fijar `--seed`.
- Registrar siempre parámetros y versión del software.
