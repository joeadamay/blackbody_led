# How to Use

## `temp_to_rgb.py`

This program produces a table of RGB colors corresponding to the emission of a
black-body at a range of temperatures and voltages.

In a terminal, run (exclude the "$") `$ python3 temp_led.py`.

The program will ask you to select a mode for the independent variable.  For
example, inputting "Voltage" or "V" will allow you to control the range of
voltages that the program uses when computing the colors.

The program will then ask for a minimum and maximum temperature as well as a
step size.  Please provide it with reasonable values (i.e. Positive values and a
minimum that is less than the maximum).  Also, try to keep the minimum values
somewhat high (i.e. a temperature above 273 K and a voltage above 10 V); the
math tends to fail when the values are too small.

### Requirements

- Python 3, with the packages:
  - csv
  - math
  - numpy
  - matplotlib (optional; this is just for visualization)
- The CIE-XYZ 1964 visible spectrum datatable, named `CIE_xyz_1964_10deg.csv`
  and located in the same directory (folder).  It can be downloaded from the
  [CIE Website](https://cie.co.at/datatable/cie-1964-colour-matching-functions-10-degree-observer "CIE 1964 Colour-Matching Functions: 10 Degree Observer").

### Output

If the program runs to completion, it will ask that you provide a destination
for the output file.  The extension `.csv` is automatically appended.

The output file is is a table in CSV format comparing absolute temperature
to the corresponding RGB color (and includes voltages if Voltage Mode was
selected).

## `blackbody_color_demo.html`

This program uses a linear approximation to display the colors of the black-body
spectrum.  To run, open the file in a browser (e.g. Chrome, Firefox, or Safari).

