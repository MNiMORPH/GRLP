# References

## GRLP papers

The model, its correction, and its published extensions and applications. All are
open access (CC-BY). Source-grounded notes on each — governing equations,
parameters, and how they map to the code — are in `docs/literature/README.md` in
the repository.

- **Wickert, A. D. and Schildgen, T. F. (2019).** Long-profile evolution of
  transport-limited gravel-bed rivers. *Earth Surface Dynamics*, **7**, 17–43.
  [doi:10.5194/esurf-7-17-2019](https://doi.org/10.5194/esurf-7-17-2019).
  *The foundational paper — defines the model.*
- **Wickert, A. D. and Schildgen, T. F. (2020).** Corrigendum to "Long-profile
  evolution of transport-limited gravel-bed rivers." *Earth Surface Dynamics*,
  **7**, 17–43.
  [doi:10.5194/esurf-7-17-2019-corrigendum](https://doi.org/10.5194/esurf-7-17-2019-corrigendum).
  *Removes the valley-width-derivative term; corrected code is GRLP v1.4.1.*
- **McNab, F., Schildgen, T. F., Turowski, J. M., and Wickert, A. D. (2023).**
  Diverse responses of alluvial rivers to periodic environmental change.
  *Geophysical Research Letters*, **50**, e2023GL103075.
  [doi:10.1029/2023GL103075](https://doi.org/10.1029/2023GL103075).
  *Gain/lag linearization and equilibration timescale.*
- **McNab, F., Schildgen, T. F., Turowski, J. M., and Wickert, A. D. (2025).**
  Influence of network geometry on long-term morphodynamics of alluvial rivers.
  *Earth Surface Dynamics*, **13**, 1059–1092.
  [doi:10.5194/esurf-13-1059-2025](https://doi.org/10.5194/esurf-13-1059-2025).
  *Extends GRLP to drainage networks.*
- **Ruby, A., McNab, F., Schildgen, T. F., Wickert, A. D., and Fernandes, V. M.
  (2026).** How sediment supply, sea-level, and glacial isostatic oscillations
  drive alluvial river long-profile evolution and terrace formation.
  *AGU Advances*, **7**, e2025AV002035.
  [doi:10.1029/2025AV002035](https://doi.org/10.1029/2025AV002035).
  *Río Santa Cruz application; shelf-coupled base-level forcing.*

## Selected supporting literature

Cited by the model derivation for its physics and closures.

- **Parker, G. (1978).** Self-formed straight rivers with equilibrium banks and
  mobile bed. Part 2. The gravel river. *Journal of Fluid Mechanics*, **89**,
  127–146. *Equilibrium-width (near-threshold) channel.*
- **Wong, M. and Parker, G. (2006).** Reanalysis and correction of bed-load
  relation of Meyer-Peter and Müller using their own database.
  *Journal of Hydraulic Engineering*, **132**, 1159–1168.
  *Bed-load transport coefficient.*
- **Meyer-Peter, E. and Müller, R. (1948).** Formulas for bed-load transport.
  *Proc. 2nd Meeting IAHR*, 39–64.
- **Paola, C., Heller, P. L., and Angevine, C. L. (1992).** The large-scale
  dynamics of grain-size variation in alluvial basins, 1: Theory.
  *Basin Research*, **4**, 73–90. *Diffusive basin-filling; equilibration time.*
- **Blom, A., Arkesteijn, L., Chavarrías, V., and Viparelli, E. (2017).** The
  equilibrium alluvial river under variable flow and its channel-forming
  discharge. *JGR Earth Surface*, **122**, 1924–1948.
- **Shreve, R. L. (1966).** Statistical law of stream numbers.
  *Journal of Geology*, **74**, 17–37; and **Shreve, R. L. (1974).** Variation of
  mainstream length with basin area in river networks. *Water Resources
  Research*, **10**, 1167–1177. *Random network topology.*
- **Hack, J. T. (1957).** Studies of longitudinal stream profiles in Virginia and
  Maryland. *USGS Professional Paper* 294-B.
- **Kennett, J. P. (1982).** *Marine Geology.* Prentice-Hall. *Continental-shelf
  gradient used for shelf-coupled base level.*
