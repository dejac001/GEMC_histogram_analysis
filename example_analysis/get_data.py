import json
from histogram.writer import histogram_bins_and_values, gaussian_function


def get_data(path):
    with open(path + 'output-long_gaussian_results.json', 'r') as f:
        gaussian_results = json.load(f)
    with open(path + 'output-long_histograms.json', 'r') as f:
        histogram_results = json.load(f)
    return gaussian_results, histogram_results


def main(g_data, h_data, output_path):
    """

    :param g_data:
    :param h_data:
    :param output_path:
    :return:
    """
    histogram_bins_and_values(output_path, h_data)
    gaussian_function(output_path, g_data)


def high_temperature():
    path = '../example_data/pentanal_545K/par1/'
    g, h = get_data(path)
    main(g, h, 'pentanal_545K_par1/')


def low_temperature():
    path = '../example_data/pentanal_510K/par1/'
    g, h = get_data(path)
    main(g, h, 'pentanal_510K_par1/')


if __name__ == '__main__':
    high_temperature()
    low_temperature()
