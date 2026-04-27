# Manual de usuario — padic-grn-core

## 1. Requisitos

- Python 3.10+
- Archivo CSV de transiciones

## 2. Instalación local (modo desarrollo)

```bash
cd software-suite/padic-grn-core
pip install -e .
```

## 3. Uso básico: un ordenamiento

```bash
pgcore-score-ordering transitions.csv 13 2 "12 4 3 8 1 5 7 2 10 11 9 6 0" --output-json scores.json
```

Resultado:

- genera `scores.json` con `mu_E`, `mu_A`, `mu_I`.

## 4. Uso por lote (batch)

`orderings.csv` debe tener una columna `ordering`.

```bash
pgcore-batch-scores transitions.csv 13 2 orderings.csv --output-csv batch_scores.csv
```

Resultado:

- `batch_scores.csv` con columnas:
  - `ordering`
  - `mu_E`
  - `mu_A`
  - `mu_I`

## 5. Buenas prácticas

- usar rutas relativas limpias y evitar espacios en nombres de archivos de trabajo.
- guardar siempre los parámetros de corrida junto al resultado.
- versionar los archivos de entrada críticos.

## 6. Solución de problemas

- Error “ordering inválido”:
  - verifica que contenga todos los índices `0..N-1` sin repetición.
- Error de CSV:
  - verifica columnas `config_from,config_to`.
- Error de módulo legado no encontrado:
  - ejecutar desde el workspace donde existe el árbol fuente oficial.
