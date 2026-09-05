# Results

Nothing here is written by hand. Every file is produced by a script in `analysis/`, reading only the inputs in `data/`, so deleting this directory and running the seven scripts again rebuilds it exactly.

| file | written by | what it holds |
| --- | --- | --- |
| `figures/systematics_budget.png` | `a05_figure.py` | the four-panel summary used in the README and the result note |
| `tables/a01_library_spread.csv` | `a01_signal_versus_library_spread.py` | the attributed offset next to the spread across the published non-magnetic libraries, per passband |
| `tables/a02_parameter_systematics.csv` | `a02_stellar_parameter_systematics.py` | the derivative of h1 and h2 with respect to each stellar parameter, and the scale error each would need to produce the offset |
| `tables/a03_facular_chord.csv` | `a03_facular_chord_test.py` | the facular filling factor the transit chord would need, at each measured contrast |
| `tables/a06_h2_reliability.csv` | `a06_h2_reliability.py` | the measured offsets against the analysis systematic Maxted (2023) reports for the same quantity |
| `tables/a07_near_solar_assumption.csv` | `a07_solar_residual.py` | how far non-magnetic limb darkening moves between solar parameters and the sample mean, in one library and one passband |
| `tables/a07_solar_passband_spread.csv` | `a07_solar_residual.py` | model h1 and h2 at solar parameters in four passbands, showing how much the passband choice is worth |

The last two columns of `a07_solar_passband_spread.csv` are differences between a measurement at 579.88 nm and a model in a broad band. They are not residuals, and the file header says so. They are printed to show the size of the passband problem, which is the subject of part one of `a07`.
