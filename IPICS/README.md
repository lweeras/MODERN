### This directory contains the data files used for testing the newly proposed Interior Point, Inner Corner Score for measuring reserve compactness.

All data files contain four sections of data, each separated by a blank line:

1) Species Representation -- Each line contains three values:
    1) Integer indicating which species data relates to
    2) Integer indicating row position of cell containing species
    3) Integer indicating column position of cell containing species
2) Cost -- Each line contains three value:
    1) Integer indicating row position of cell
    2) Integer indicating column position of cell
    3) Decimal value indicating cost to select the given cell as part of the reserve
3) Required Representation -- One line, containing one value for each species. The $i^{\text{th}}$ value indicates the number of cells containing the $i^{\text{th}}$ species that must be included in the reserve.
4) Maximum Budget -- A single value indicating the maximum budget that can be spent on selecting a reserve.

Note that some directories contain a "modified" file that contains a larger budget so as to ensure that the associated instance is feasible.
