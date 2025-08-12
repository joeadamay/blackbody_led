import csv
import math
import matplotlib.pyplot as plt
import numpy as np


# Compute the spectral radiance at a given wavelength (in meters) and
# temperature (in Kelvin)
# See https://en.wikipedia.org/wiki/Planck's_law for more information.
def planck(wavelength, temperature):
    # Planck Constant (J Hz^-1)
    h = 6.626_070_15e-34
    # Speed of Light (m s^-1)
    c = 2.997_924_58e8
    # Boltzman Constant (J K)
    k_b = 1.380_649e-23

    frequency = c / wavelength

    numerator = 2.0 * h * pow(frequency, 3)
    denominator = c * c * math.expm1(h * frequency / (k_b * temperature))

    return numerator / denominator


# Use Composite Simpson's Rule to integrate over the points
# `values` is a list of 'y-values' and will be shortened if not of even length
# `subint_size` is the difference between the 'x-values' of the points. It is
#    assumed to be constant.
# See Sauer, T. (2018). *Numerical Analysis* (3rd ed.). Pearson. for more info.
def simpson(values, subint_size):
    num_points = len(values)
    # Check that values has an even length
    if num_points % 2 != 0:
        num_points -= 1

    num_panels = int(num_points / 2)
    result = values[0] + values[2 * num_panels]

    sum_value = 0.0
    for i in range(1, num_panels + 1):
        sum_value += values[2 * i - 1]
    result = 4.0 * sum_value

    sum_value = 0.0
    for i in range(1, num_panels):
        sum_value += values[2 * i]
    result += 2.0 * sum_value

    return result


# `temperature` is in Kelvin
# `data` is a list of lists describing the perceptive response of the visible
#  spectrum: [wavelength (nm), CIE X, CIE Y, CIE Z]
def get_rgb_from_temp(temperature, data):
    # Get the black-body spectrum
    blackbody = []
    for datum in data:
        wavelength = datum[0] * 1e-9 # Convert from nm to meters
        radiance = planck(wavelength, temperature)
        blackbody.append(radiance)

    # Compute the product between XYZ colorspace and black-body radiation
    x_bar = [data[i][1] * blackbody[i] for i in range(len(data))]
    y_bar = [data[i][2] * blackbody[i] for i in range(len(data))]
    z_bar = [data[i][3] * blackbody[i] for i in range(len(data))]

    # Compute the integral across the region of spectrum represented by the data
    # points
    # Get the subinterval size, assuming that all points are evenly spaced
    subint_size = data[1][0] - data[0][0]
    x = simpson(x_bar, subint_size)
    y = simpson(y_bar, subint_size)
    z = simpson(z_bar, subint_size)
    xyz = np.array([x, y, z])

    # I don't remember where I got this from: probably from Wikipedia.
    # Let's hope it works.
    xyz_to_rgb = np.array([[ 1.91020, -1.11212,  0.20191],
                           [ 0.37095,  0.62905,  0.0    ],
                           [ 0.0    ,  0.0    ,  1.0    ]])
    xyz_to_rgb = np.linalg.inv(xyz_to_rgb)

    # Convert from the CIE-XYZ color representation to RGB
    rgb = xyz_to_rgb @ xyz
    # Normalize the color to have maximal lightness (i.e. the vector is scaled
    # so that the greatest component becomes 1.0)
    rgb_normalized = rgb / np.linalg.norm(rgb, np.inf)
    #print(f'Temperature: {temperature}')
    #print(f'\tXYZ: {xyz}')
    #print(f'\tRGB: {rgb}')
    #print(f'\tRGB Normalized: {rgb_normalized}')

    return rgb_normalized


if __name__ == '__main__':
    input_file_name = 'CIE_xyz_1964_10deg.csv'
    #min_temp =  1_000
    #max_temp = 10_000
    #temp_step =   100
    try:
        min_temp =  float(input('Minimum temperature (K): '))
        max_temp =  float(input('Maximum temperature (K): '))
        temp_step = float(input('Step size (K): '))

        assert(min_temp > 0)
        assert(max_temp > 0)
        assert(temp_step > 0)
        assert(max_temp > min_temp)
    except:
        print("\nPlease provide decimal numbers with reasonable values,")
        print("i.e. positive temperatures and a minimum that is less than the")
        print("maximum.\n")
        assert(False)

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

    # Get Temperature-RGB relationship
    temp_rgb = []
    temp = min_temp
    while temp <= max_temp:
        rgb = get_rgb_from_temp(temp, data)
        temp_rgb.append([temp, rgb])
        temp += temp_step

    # Write to file
    output_file_name = input('Output file name: ')
    output_file_name += ".csv"
    with open(output_file_name, 'w') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(['Temperature (K)', 'Red', 'Green', 'Blue'])
        for row in temp_rgb:
            temp = row[0]
            r = row[1][0]
            g = row[1][1]
            b = row[1][2]
            writer.writerow([temp, r, g, b])

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

