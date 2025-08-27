import csv
import math
import matplotlib.pyplot as plt
import numpy as np


# Compute the spectral radiance (W sr^-1 m^-3) at a given wavelength (in meters)
# and temperature (in Kelvin)
# For more information, see https://en.wikipedia.org/wiki/Planck's_law or the
# [CIE Technical Report Colorimetry](https://archive.org/details/gov.law.cie.15.2004)
def planck(wavelength, temperature):
    # Planck Constant (J Hz^-1)
    h = 6.626_070_15e-34
    # Speed of Light (m s^-1)
    c = 2.997_924_58e8
    # Boltzman Constant (J K^-1)
    k_b = 1.380_649e-23

    numerator = 2.0 * h * c * c
    denominator = pow(wavelength, 5) * math.expm1(h * c / (wavelength * k_b * temperature))

    return numerator / denominator


# Compute the temperature of a incandescent lamp's filament, given a voltage
# (in Volts).  The empirical method used here is due to Martin Kykta.
# Note that this relation only holds if the voltage is not near zero.
# See https://pubs.aip.org/aip/adv/article/12/10/105116/2819829/Incandescent-lamp-design-and-lifetime
# for more information.
def voltage_to_temp(voltage):
    if voltage < 0:
        voltage *= -1

    # Stefan-Boltzmann Constant (from
    # https://en.wikipedia.org/wiki/Stefan%E2%80%93Boltzmann_law)
    # In W m^-2 K^-4
    stef_boltz = 5.670_374_419e-8
    # Emissivity of tungsten (from Kykta, 2022)
    emissivity = 0.28
    # In meters
    #filament_length = 0.02814
    filament_length = 0.023_14 # GE47 from Kykta, 2022
    # In meters
    #filament_radius = 11.5533e-6
    filament_radius = 0.000_010_91 # GE47 from Kykta, 2022

    b_1 = (2.0 * stef_boltz * emissivity * 2.96e8 ** 4) ** -0.232
    b_1 /= math.pi

    b_2 = (2.0 * math.pi * stef_boltz * emissivity * b_1) ** -0.25

    B_2 = b_2 * filament_length ** -0.384 * filament_radius ** 0.192

    temperature = B_2 * voltage ** 0.384

    return temperature


# Use Composite Simpson's Rule to integrate over the points
# `values` is a list of 'y-values' and will be shortened if not of even length
# `subint_size` is the difference between the 'x-values' of the points. It is
#    assumed to be constant.
# See Sauer, T. (2018). *Numerical Analysis* (3rd ed.). Pearson. for more info.
def simpson(values, subint_size):
    num_points = len(values)
    # Check that there are an even number of subintervals
    if (num_points - 1) % 2 != 0:
        num_points -= 1

    num_panels = int(num_points / 2)

    result = values[0] + values[2 * num_panels]

    sum_value = 0.0
    for i in range(1, num_panels + 1):
        sum_value += values[2 * i - 1]
    result += 4.0 * sum_value

    sum_value = 0.0
    for i in range(1, num_panels):
        sum_value += values[2 * i]
    result += 2.0 * sum_value

    result *= subint_size / 3.0

    return result


# `temperature` is in Kelvin
# `data` is a list of lists describing the perceptive response of the visible
#   spectrum: [wavelength (nm), CIE X, CIE Y, CIE Z]
# Returns xyz
def get_xyz_from_temp(temperature, data):
    x_bar = []
    y_bar = []
    z_bar = []
    for i in range(len(data)):
        wavelength = data[i][0] * 1e-9 # Convert from nm to meters
        # Get the blackbody spectral radiance
        spectral_radiance = planck(wavelength, temperature)
        # Multiply the XYZ color-matching functions by the spectral radiance
        x_bar.append(data[i][1] * spectral_radiance)
        y_bar.append(data[i][2] * spectral_radiance)
        z_bar.append(data[i][3] * spectral_radiance)

    # Compute the integral across the region of spectrum represented by the data
    # points
    # Get the subinterval size, assuming that all points are evenly spaced
    subint_size = data[1][0] - data[0][0]
    x = simpson(x_bar, subint_size)
    y = simpson(y_bar, subint_size)
    z = simpson(z_bar, subint_size)
    xyz = np.array([x, y, z])

    # Multiply by maximum spectral luminous efficacy, as per CIE guidelines, to
    # go from radiance (W sr^-1 m^-2) to luminance (lm sr^-1 m^-2 or cd m^-2)
    xyz *= 683 # lm W^-1

    return xyz


# `xyz` is a numpy array
def xyz_to_rgb(xyz):
    # From CIE Colorimetry 3rd Edition (2004)
    # https://archive.org/details/gov.law.cie.15.2004
    xyz_to_rgb = np.array([[2.768_892, 1.751_748, 1.130_160],
                           [1.000_000, 4.590_700, 0.060_100],
                           [0.0      , 0.056_508, 5.594_292]])
    xyz_to_rgb = np.linalg.inv(xyz_to_rgb)

    # Convert from the CIE-XYZ color representation to RGB
    rgb = xyz_to_rgb @ xyz

    return rgb


def main():
    input_file_name = 'CIE_xyz_1964_10deg.csv'

    try:
        selection = input('Mode ((T)emperature/(V)oltage: ').lower()
        voltage_mode = False
        if selection[0] == 't':
            voltage_mode = False
        elif selection[0] == 'v':
            voltage_mode = True
        else:
            assert(False)
    except:
        print('\nPlease answer with "T" or "V".\n')
        return

    try:
        value_text = 'voltage (V)' if voltage_mode else 'temperature (K)'
        min_value           =  float(input(f'Minimum {value_text}: '))
        max_value           =  float(input(f'Maximum {value_text}: '))
        step_size           = float(input('Step size: '))
        reference_x         = float(input(f'Reference {value_text}: '))
        reference_luminance = float(input(f'Reference Luminance (lm sr^-1 m^-2): '))

        assert(min_value > 0)
        assert(max_value > 0)
        assert(step_size > 0)
        assert(max_value > min_value)
        assert(reference_x > 0)
        assert(reference_luminance > 0)
    except:
        print("\nPlease provide decimal numbers with reasonable values,")
        print("i.e. positive and a minimum that is less than the maximum.\n")
        return

    # Start reading the file
    in_f = open(input_file_name, 'r')
    data_reader = csv.reader(in_f)

    data = []

    # Extract and clean data
    for datum in data_reader:
        # Remove bad data
        if len(datum) != 4:
            print(f'Excluding bad datum: {datum}')
            continue

        # Convert strings to numbers. Assume nonnumerics (e.g. 'NaN') are zero.
        for j in range(len(datum)):
            try:
                datum[j] = float(datum[j])
                # Don't include +/-Inf and NaN
                assert(math.isfinite(datum[j]))
            except:
                datum[j] = 0.0

        data.append(datum)

    # Close the file
    in_f.close()

    # Calculate the ratio between the true luminance and the calculated
    # luminance at the reference point (The Y of XYZ is luminance)
    ref_temp = voltage_to_temp(reference_x) if voltage_mode else reference_x
    ref_xyz = get_xyz_from_temp(ref_temp, data)
    calibration_coefficient = reference_luminance / ref_xyz[1]

    # Get Temperature-RGB relationship
    voltages = []
    temp_xyz_rgb = []
    cur_value = min_value
    while cur_value <= max_value:
        # Calculate the temperature if necessary
        temp = voltage_to_temp(cur_value) if voltage_mode else cur_value

        xyz = get_xyz_from_temp(temp, data)
        xyz *= calibration_coefficient
        rgb = xyz_to_rgb(xyz)
        temp_xyz_rgb.append([temp, xyz, rgb])

        # Record the voltage if necessary
        if voltage_mode:
            voltages.append(cur_value)

        cur_value += step_size

    # Write to file
    output_file_name = input('Output file name: ')
    output_file_name += ".csv"
    with open(output_file_name, 'w') as out_f:
        header = []
        if voltage_mode:
            header.append('Voltage (V)')
        header.extend(['Temperature (K)', 'X', 'Y', 'Z', 'Red', 'Green', 'Blue'])

        # Add the reference point
        reference_csv = []
        if voltage_mode:
            reference_csv.append('Reference Voltage (V):')
        else:
            reference_csv.append('Reference Temperature (K):')
        reference_csv.append(reference_x)
        reference_csv.append('Reference Luminance (lm sr^-1 m^-2):')
        reference_csv.append(reference_luminance)

        writer = csv.writer(out_f)
        writer.writerow(reference_csv)
        writer.writerow(header)
        for i in range(len(temp_xyz_rgb)):
            row = temp_xyz_rgb[i]

            temp = row[0]
            [x, y, z] = row[1]
            [r, g, b] = row[2]

            data = []
            # Include voltage if necessary
            if voltage_mode:
                data.append(voltages[i])
            # Include the temperature and color
            data.extend([temp, x, y, z, r, g, b])
            writer.writerow(data)

    # Make a plot
    #x_axis = [d[0] for d in temp_rgb]
    #r = [d[1][0] for d in temp_rgb]
    #g = [d[1][1] for d in temp_rgb]
    #b = [d[1][2] for d in temp_rgb]
    #rgb = [[r[i], g[i], b[i]] for i in range(len(r))]

    #plt.scatter(x_axis, r, c=rgb, marker='.', label='Red')
    #plt.scatter(x_axis, g, c=rgb, marker=',', label='Green')
    #plt.scatter(x_axis, b, c=rgb, marker='o', label='Blue')

    ##plt.plot(x_axis, r, color='red', label='Red')
    ##plt.plot(x_axis, g, color='green', label='Green')
    ##plt.plot(x_axis, b, color='blue', label='Blue')

    #plt.xlabel('Temperature (K)')
    #plt.ylabel('Standardized Color Value')
    #plt.legend()
    #plt.gca().set_facecolor('black')
    #plt.show()


if __name__ == '__main__':
    main()

