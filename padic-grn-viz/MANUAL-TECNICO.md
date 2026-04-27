# Manual técnico — padic-grn-viz

## 1. Propósito

`padic-grn-viz` transforma resultados numéricos (`core/search`) en artefactos visuales reproducibles para:

- análisis,
- suplementario de artículos,
- demostración web institucional.

## 2. Arquitectura actual

- `src/padic_grn_viz/engine.py`
  - carga de resultados (`load_scores_json`),
  - gráfico de barras de métricas (`plot_score_bars`),
  - gráfico de ordering (`plot_ordering_positions`).
- `src/padic_grn_viz/cli.py`
  - `pgviz-dependency-graph`
  - `pgviz-hierarchical-graph`

## 3. Tipos de salida

- figuras `png/pdf`,
- resúmenes tabulares de apoyo.

## 4. Contrato mínimo recomendado para integración

- entradas en `csv/json` con:
  - metadata (`N`, `p`, red, fecha),
  - métricas (`mu_E`, `mu_A`, `mu_I`),
  - ordering usado.

## 5. Roadmap técnico

1. Añadir gráficos avanzados de perfiles por nivel (`E/A/I`).
2. Añadir plantillas de reporte visual unificado (PNG + PDF + resumen).
3. Integración directa con backend de demo institucional.
