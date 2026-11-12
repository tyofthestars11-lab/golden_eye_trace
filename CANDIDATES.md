# Candidate table — raw and corrected

Every candidate from `prereg/ovro_survey.json`, with the corrections that can be applied without a sky position. `earth max` is the largest drift Earth's rotation and orbit can induce at that frequency, over the whole sky, for that date; `dv/dt` is the acceleration a transmitter would need to produce the measured drift; `off` is whether an off-source observation exists in this product (it does not).

| date | MHz | drift Hz/s | SNR | res Hz/s | earth max Hz/s | dv/dt m/s² | band | repeats | off | k | lock dist | φ phase° | disposition |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 20240207 | 54.4530 | 8,194.7 | 7.68 | 245.5 | 0.01 | 45,116 | inside the 54-60 MHz TV channel 2 allocation | 0 | none | 1 | 0.945 | 110.5 | UNRESOLVED (in a TV band, but no repeat across days) |
| 20240207 | 53.7833 | 6,970.9 | 7.16 | 392.8 | 0.01 | 38,856 | inside the 50-54 MHz amateur 6 m allocation | 0 | none | 1 | 0.954 | 101.2 | UNRESOLVED |
| 20240207 | 55.9010 | 7,866.9 | 7.00 | 245.5 | 0.01 | 42,189 | inside the 54-60 MHz TV channel 2 allocation | 0 | none | 1 | 0.948 | 130.1 | UNRESOLVED (in a TV band; repeat untestable, outside survey coverage) |
| 20240322 | 47.2416 | -172.8 | 6.67 | 392.8 | 0.00 | 1,097 | below the 6 m band; no allocated broadcast service here | 1 | none | 1 | 0.999 | 4.2 | UNRESOLVED |
| 20240324 | 47.3143 | -504.0 | 11.67 | 392.8 | 0.00 | 3,193 | below the 6 m band; no allocated broadcast service here | 1 | none | 1 | 0.997 | 5.3 | UNRESOLVED |
| 20240325 | 45.9843 | 41,190.2 | 953.98 | 392.8 | 0.00 | 268,538 | below the 6 m band; no allocated broadcast service here | 0 | none | 1 | 0.725 | 344.0 | ASTROPHYSICAL (bright burst group) |
| 20240325 | 46.9828 | 41,185.4 | 604.03 | 392.8 | 0.00 | 262,800 | below the 6 m band; no allocated broadcast service here | 0 | none | 1 | 0.725 | 0.1 | ASTROPHYSICAL (bright burst group) |
| 20240325 | 47.9380 | 41,184.1 | 447.08 | 392.8 | 0.00 | 257,555 | below the 6 m band; no allocated broadcast service here | 0 | none | 1 | 0.725 | 15.1 | ASTROPHYSICAL (bright burst group) |
| 20240325 | 48.9937 | 41,179.1 | 171.23 | 392.8 | 0.00 | 251,975 | below the 6 m band; no allocated broadcast service here | 0 | none | 1 | 0.726 | 31.4 | ASTROPHYSICAL (bright burst group) |
| 20240920 | 45.4206 | 195,015.7 | 8.66 | 392.8 | 0.00 | 1,287,173 | below the 6 m band; no allocated broadcast service here | 0 | none | 1 | 0.300 | 334.8 | INSTRUMENTAL ALIASING (rate at Nyquist) |
| 20240920 | 46.5890 | 194,896.8 | 7.91 | 392.8 | 0.00 | 1,254,130 | below the 6 m band; no allocated broadcast service here | 0 | none | 1 | 0.299 | 353.8 | INSTRUMENTAL ALIASING (rate at Nyquist) |
| 20240920 | 47.5358 | 194,905.6 | 6.49 | 392.8 | 0.00 | 1,229,206 | below the 6 m band; no allocated broadcast service here | 0 | none | 1 | 0.299 | 8.8 | INSTRUMENTAL ALIASING (rate at Nyquist) |

Phi control: best lock distance 0.2991 against 2,000 randomized-phase draws, p = 1.0000. The phi columns are a preregistered feature, not evidence.
