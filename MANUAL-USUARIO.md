# Manual de usuario — p-adic GRN Suite

## 1. Qué es p-adic GRN Suite

`p-adic GRN Suite` integra tres módulos:

- `padic-grn-core`: cálculo de métricas p-ádicas (`mu_E`, `mu_A`, `mu_I`).
- `padic-grn-search`: búsqueda/certificación de orderings.
- `padic-grn-viz`: visualización y reporte gráfico.

Flujo recomendado de uso:

1. calcular métricas con `core`,
2. explorar/certificar orderings con `search`,
3. visualizar y documentar resultados con `viz`.

## 2. Requisitos

- Python 3.10+
- `pip`
- entorno local con acceso a los tres submódulos en `software-suite/`

## 3. Instalación (modo desarrollo)

```bash
cd software-suite/padic-grn-core && pip install -e .
cd ../padic-grn-search && pip install -e .
cd ../padic-grn-viz && pip install -e .
```

## 4. Entrada mínima esperada

Para `core/search` se usa un CSV de transiciones discretas con columnas:

- `config_from`
- `config_to`

Para `viz` se usan resultados JSON/CSV producidos por las corridas anteriores.

## 5. Uso rápido de la suite

### 5.1 Calcular un ordering con core

```bash
cd software-suite/padic-grn-core
pgcore-score-ordering examples/transitions_p2_n2_identity.csv 2 2 "0 1" --output-json core_scores.json
```

Salida esperada: archivo `core_scores.json`.

### 5.2 Búsqueda con search

```bash
cd ../padic-grn-search
pgsearch-ga ../padic-grn-core/examples/transitions_p2_n2_identity.csv 2 2 --objective minimize_mu_E --generations 5 --population 10
```

Salida esperada: mejor ordering y métricas objetivo.

### 5.3 Visualización con viz

```bash
cd ../padic-grn-viz
pgviz-dependency-graph
pgviz-hierarchical-graph
```

Salida esperada: archivos de figura/reportes (según parámetros de ejecución).

## 6. Ejecución de validación integrada

```bash
cd software-suite
./RUN-ALL-VALIDATIONS.sh
```

Salida esperada:

- pruebas `pytest` para `core`, `search`, `viz`,
- estado global de validación de la suite.

## 7. Buenas prácticas

- guardar junto al resultado: parámetros, versión y fecha de corrida,
- mantener nombres de archivos estables para reproducibilidad,
- usar siempre ejemplos pequeños antes de correr casos grandes.

## 8. Solución de problemas

- **Comando no encontrado**: verificar instalación editable (`pip install -e .`) en cada módulo.
- **CSV inválido**: revisar columnas y formato de tabla de transiciones.
- **Ordering inválido**: verificar que sea una permutación sin repetición de `0..N-1`.
- **Error de visualización**: confirmar instalación de dependencias gráficas.
