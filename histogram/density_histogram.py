"""
"""

from __future__ import division, print_function

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from .chem_constants import N_av, R
from . import reader
import json


def gaussian(x, a, b, c):
    """Gaussian function

    .. math::
        f(x) = a e^{-(x-b)^2 / 2 / c^2}

    :param x: input variable
    :param a: height of curve peak
    :param b: center of peak
    :param c: standard deviation
    :return: f(x)
    """
    square_term = (x - b) / c
    exponential_term = np.exp(-1. / 2. * square_term * square_term)

    return a * exponential_term


def jac_gaussian(x, a, b, c):
    """Jacobian of :py:func:`~histogram.density_histogram.gaussian`

    :param x: input variable
    :param a: height of curve peak
    :param b: center of peak
    :param c: standard deviation
    :return:
        .. math:: \\frac{\partial f}{\partial a}, \\frac{\partial f}{\partial b},
            \\frac{\partial f}{\partial c}
    """
    ds = (x - b)
    square_term = ds / c
    power_term = -1. / 2. * square_term * square_term
    exponential_term = np.exp(power_term)
    f = a * exponential_term

    J = np.array([
        exponential_term,
        f * ds / c / c,
        f * ds * ds / c / c / c
    ])

    return J.transpose()


def calculate_density(volumes, number_of_chains, cycles):
    """Convert input number of chains and volumes to densities

    .. warning:
        Will be incorrect for boxes where molecules have biasing potentials.
        Does not correct the number density for biasing potential in box.

    .. note:
        Biasing potential could be compatible by reading values in run1a.dat file and storing
        in :py:class:`~histogram.density_histogram.Blocks`

    :param volumes: volume of box at end of each cycle in nm**3
    :param number_of_chains: number of molecules of each type and in each box during each cycle
    :type number_of_chains: dict
    :param cycles: cycles over which to calculate density
    :type cycles: list
    :type volumes: dict
    :return: densities
    """
    density = {}
    for mol in number_of_chains.keys():
        density[mol] = {}
        for box in number_of_chains[mol].keys():
            density[mol][box] = {}
            for cycle in cycles:
                assert cycle in volumes[box].keys(), 'Cycle %i not in Volumes' % cycle
                assert cycle in number_of_chains[mol][box].keys(), 'Cycle %i not in number of chains' % cycle
                density[mol][box][cycle] = number_of_chains[mol][box][cycle] / volumes[box][cycle]
    return density


def calculate_compressibility(total_densities, pressures, T):
    """Calculate compressibility factors from total densities and pressures

    :param total_densities: total densities of boxes in molecules / nm**3
    :type total_densities: dict
    :param pressures: pressures in kPa
    :type pressures: dict
    :param T: Temperature in Kelvin
    :type T: float
    :return: List of compressibility factors of each box for each point provided
    """
    from .chem_constants import N_av, R
    Z = {}
    for box in pressures.keys():
        Z[box] = {}
        for cycle in pressures[box].keys():
            Z[box][cycle] = pressures[box][cycle]/total_densities[box][cycle] * N_av / R['nm**3*kPa/(mol*K)']/T
    return Z


def make_block_averages(data, number_blocks):
    """Make block averages based off of number of blocks and data provided.

    :param data: a thermodynamic property measured at each point for each cycle
    :type data: dict
    :param number_blocks: number of nodes over which to average property
    :type number_blocks: int
    :return: means and standard deviations of data for each box
    """
    block_averages = {}
    for box in data.keys():
        vals = list(data[box].values())
        block_averages[box] = []
        number_cycles = len(vals)
        cycles_per_block = number_cycles // number_blocks
        for i in range(number_blocks):
            my_values = vals[i * cycles_per_block:(i * cycles_per_block + cycles_per_block)]
            block_averages[box].append(np.mean(my_values))
    data = {}
    for box, vals in block_averages.items():
        data[box] = {'mean': np.mean(vals), 'stdev': np.std(vals)}
    return data


def calculate_total_densities(dens, boxes):
    """Calculate the total (molar) density in a box by summing individual number densities

    .. math:
        \rho_{\mathrm{total}} = \sum_i\rho_i

    :param dens: densities in molecule/nm**3 for each box for each molecule
    :type dens: dict
    :param boxes: boxes of simulation
    :type boxes: list
    :return:
    """
    total_densities = {}
    for box in boxes:
        total_densities[box] = {}
        for cycle in dens['molecule 1'][box].keys():
            my_total_dens = 0.
            for molty in dens.keys():
                my_total_dens += dens[molty][box][cycle]
            total_densities[box][cycle] = my_total_dens
    return total_densities


def update_data_box(store_data, new_data):
    """Updata information stored for data that is only box specific

    :param store_data: stored information of given data
    :param new_data: new information to add to stored information
    :return: None
    """
    for box in new_data.keys():
        if box not in store_data.keys():
            store_data[box] = {}
        for cycle, val in new_data[box].items():
            assert cycle not in store_data[box].keys(), 'Cycle {} already in stored data'.format(cycle)
            store_data[box][cycle] = val


def update_data_molty(store_data, new_data):
    """Updata information stored for data that depends on both box and molty

    :param store_data: stored information of given data
    :param new_data: new information to add to stored information
    :return: None
    """
    for nmolty in new_data.keys():
        if nmolty not in store_data.keys():
            store_data[nmolty] = {}
        update_data_box(store_data[nmolty], new_data[nmolty])


class Blocks:
    """
    Read input arguments and initialize data types
        for one independent simulation.
    """

    def __init__(self, **kwargs):

        self.number_of_blocks = 10  # only matters for computing error of ONE independent simulation
        self.files = kwargs['files']
        self.T = kwargs['temperature']
        self.unit = kwargs['d_units']
        self.output_file_name = kwargs['output_file_name']
        Blocks.check_arguments(self)
        self.units = {
            'pressure': 'kPa',
            'density': 'chains / nm**3',
            'compressibility': 'dimensionless'
        }
        self.data = None
        self.total_densities = None
        self.pressures = {}
        self.compressibilities = None
        self.block_averages = {
            'pressure': {},
            'compressibility': {},
            'density': {}
        }
        self.total_densities_with_pressure = None
        self.iratp = False
        self.output_path = None
        self.molecular_weight = None

    def choose_output_file_path(self):
        """Based off of input files, try to guess where to output

        .. warning: This feature has not been tested thoroughly on a wide distribution
            of possible input arguments. It might put output files
            in places that are not expected

        :return: None
        """
        fort_12_paths = []
        for file in self.files:
            if '/' in file:
                path = file[:file.rfind('/')] + '/'
            else:
                path = './'
            if path not in fort_12_paths:
                fort_12_paths.append(path)
        if len(fort_12_paths) == 1:
            self.output_path = fort_12_paths[0]
        else:
            print('   ambiguous output path, choosing the first one')
            self.output_path = fort_12_paths[0]

    def check_arguments(self):
        """
        Check that arguments passed in are realistic
        :return: None
        """
        assert self.files, 'At least one trajectory file needs to be provided'
        assert self.T, 'Temperature needs to be input'

    def initialize_data(self):
        """Read all trajectory files and input into class attribute data types

        :return: None
        """
        all_data = {key: {} for key in ['N', 'V']}
        total_cycles = 0
        for file in self.files:
            (m_weights, f_nmolty, f_volume,
             f_pressure, iratp, total_cycles) = reader.read_fort12(file, total_cycles)
            if not self.iratp:
                self.iratp = iratp
                self.molecular_weight = m_weights['molecule 1']  #
            update_data_molty(all_data['N'], f_nmolty)
            update_data_box(all_data['V'], f_volume)
            update_data_box(self.pressures, f_pressure)
        boxes = list(self.pressures.keys())
        all_densities = calculate_density(all_data['V'], all_data['N'],
                                          cycles=list(all_data['V']['box 1'].keys()))
        self.total_densities = calculate_total_densities(all_densities, boxes)
        densities_with_pressure = calculate_density(all_data['V'], all_data['N'],
                                                    cycles=list(self.pressures['box 1'].keys()))
        self.total_densities_with_pressure = calculate_total_densities(densities_with_pressure, boxes)
        self.compressibilities = calculate_compressibility(self.total_densities_with_pressure, self.pressures, self.T)
        self.block_averages['pressure'] = make_block_averages(self.pressures, self.number_of_blocks)
        self.block_averages['compressibility'] = make_block_averages(self.compressibilities, self.number_of_blocks)
        self.block_averages['density'] = make_block_averages(self.total_densities, self.number_of_blocks)

    def output_blocks(self):
        """
        Output block averages into a csv file or pprint based off of input arguments.

        Convert units to g/mL if this was requested in args

        :return: None
        """
        if (self.unit == "g/mL"):
            self.units["density"] = "g/mL"
            vals = self.block_averages["density"]
            vals['box 1']['mean'] *= (1 / N_av) * (1 / 1e-21) * self.molecular_weight
            vals['box 2']['mean'] *= (1 / N_av) * (1 / 1e-21) * self.molecular_weight
            vals['box 1']['stdev'] *= (1 / N_av) * (1 / 1e-21) * self.molecular_weight
            vals['box 1']['stdev'] *= (1 / N_av) * (1 / 1e-21) * self.molecular_weight
        if not self.output_file_name:
            import pprint
            pprint.pprint(self.data)
        else:
            with open(self.output_path + self.output_file_name + '_block_averages.csv', 'w') as g:
                g.write('Property,Units,Box 1 Density Mean,Box 1 Density Standard Deviation,'
                        'Box 2 Density Mean,Box 2 Density Standard Deviation\n')
                for property in 'pressure', 'compressibility', 'density':
                    vals = self.block_averages[property]
                    g.write('%s,%s,%e,%e,%e,%e\n' % (property, self.units[property],
                                                     vals['box 1']['mean'], vals['box 1']['stdev'],
                                                     vals['box 2']['mean'], vals['box 2']['stdev']
                                                     ))

    def main(self):
        """
        Main method of :class:
        First, calls :py:func:`~histogram.density_histogram.Blocks.choose_output_file_path`
        Then, calls :py:func:`~histogram.density_histogram.Blocks.initialize_data`
        Finally, outputs data from
        :py:func:`~histogram.density_histogram.Blocks.output_blocks`
        :return: None
        """
        self.choose_output_file_path()
        self.initialize_data()
        self.output_blocks()


class DPD(Blocks):
    """Read input arguments and initialize data types
    for one independent simulation.
    Then, make associated histograms, plot them if desired,
    and output them into desired format.
    """
    def __init__(self, **kwargs):
        """

        :param kwargs: input arguments

        :py:attr:`~histograms.density_histogram.region_data_above_threshold` is a nested dictionary
            that contains all final information within the histogram
            after the gaussian fits of the peak.
            The data structure is {i: {j: k: value}}}
            where i is pressure, compressibility, or density,
            j is low dens or high dens (the different region names),
            k is the string 'a-b' where a is the cycle number and b is the box type,
            and value is the value at i, j, k
        """
        Blocks.__init__(self, **kwargs)
        self.figure_number = 0  # for plotting results with matplotlib
        self.fraction_peak_maximum_density = kwargs['fraction_peak_maximum_density']
        self.number_of_bins = kwargs['number_of_bins']
        self.make_plot = kwargs['plot']
        self.unit = kwargs['d_units']
        self.region_names = ('low dens',  # lower density name
                             'high dens')  # higher density name
        self.region_data_above_threshold = {
            'pressure': {key: {} for key in self.region_names},
            'compressibility': {key: {} for key in self.region_names},
            'density': {key: {} for key in self.region_names}
        }
        self.region_data = {
            'pressure': {key: {} for key in self.region_names},
            'compressibility': {key: {} for key in self.region_names},
            'density': {key: {} for key in self.region_names}
        }
        self.histograms = {
            'pressure': {},
            'compressibility': {},
            'density': {}
        }
        self.gaussian_fit_results = {}
        self.region_membership = {}
        self.check_arguments()

    def check_arguments(self):
        """
        Check if arguments passed in are realistic
        :return: None
        """
        assert self.fraction_peak_maximum_density > 0.5, 'Why use a density fraction of peak maximum less than 0.5'

    def fit_gaussian(self, bin_means, values):
        """Fit gaussian function, :py:func:`~histogram.density_histogram.gaussian`

        :param bin_means: means of density bins (nodes)
        :param values: values at means of density bins

        .. note:
            If less than three indices are above threshold, cannot fit gaussian.
            If less than ten indices are above threshold, prints warning

        The initial guess for the peak height, mean, and standard deviation
            of the gaussian are obtained from the maximum of the values
            above the threshold,
            the mean of the densities above the threshold,
            and the standard deviation of the
            densities above the threshold, respectively.

        The minimum value allowed for the peak height, mean, and standard
            deviation of the gaussian are 0, the minimum density
            above the threshold, and 0, respectively.

        The maximum value allowed for the peak height and mean
            of the gaussian are twice the maximum
            of the values above the threshold and the maximum density
            above the threshold, respectively.
            The standard deviation of the gaussian is not
            constrained to a finite maximum threshold.

        The parameters are fit using :py:func:`curve_fit` with the
            analytical jacobian :py:func:`~histogram.density_histogram.J_gaussian`
            and the initial guess and bounds listed above.

        Finally, the bounds for the values according to the fraction of the threshold
            are converted into continuous bounds in density.

        :returns:
            - params(list) the gaussian parameters fit
            - covariance(list) is the covariance of the gaussian parameters
            - min_rho_gauss(float) is the minimum density of the gaussian within the threshold specified
            - max_rho_gauss(float) is the minimum density of the gaussian within the threshold specified
            - num_indices(int) the number of indices used in the fitting of the gaussian
        """
        assert len(bin_means) > 2, 'Not enough indices for gaussian fitting, consider increasing bins'
        if len(bin_means) < 10:
            print(' Only %i indices used in fitting gaussian ' % (len(bin_means)))
        # initial guess for a, b, c = peak_height, mean, stdev
        initial_guess = [
            np.max(values), np.mean(bin_means), np.std(list(i*j for i,j in zip(bin_means, values)))
        ]
        # bounds on each variable can help
        lower_bounds = [0., np.min(bin_means), 0.]
        upper_bounds = [2.0 * np.max(values), np.max(bin_means), np.inf]
        params, covariance = curve_fit(gaussian, bin_means, values, p0=initial_guess,
                                       check_finite=True, jac=jac_gaussian,
                                       bounds=[lower_bounds, upper_bounds], method='trf',
                                       **{'ftol': 1e-12, 'xtol': 1e-12, 'gtol': 1e-12})

        # quadratic equation to determine minimum within fraction threshold
        b = -2 * params[1]
        c = params[1] * params[1] + 2 * params[2] * params[2] * np.log(self.fraction_peak_maximum_density)
        a = 1
        min_rho_gauss = (-b - np.sqrt(b * b - 4 * a * c)) / 2. / a
        max_rho_gauss = (-b + np.sqrt(b * b - 4 * a * c)) / 2. / a
        return params.tolist(), covariance.tolist(), min_rho_gauss, max_rho_gauss

    def fit_density_peaks(self):
        """Fit gaussian to both density peaks. The following procedure is followed

        1. Convert total box densities for each box type into one single array.

        2. Calculate the mean density out of all densities observed.

        3. Partition the densities into two regions--those which
            have values lower or higher than this mean value.

        4. Within each new distribution,
            perform a histogram using the number of bins requested
            during instantiation. The edges of the histogram
            correspond to the number of bins chosen, and the values
            are their corresponding probabilities normalized
            so that the sum of all probabilities is 1.

        5. Partition histogram bins into those which have values
            above the fractional threshold of the maximum value.
            Fit a gaussian function of the probabilities of
            observing each density as a function of density.

        Store the results of the binning for each region in
        :py:attr:`~histograms.density_histogram.DPD.histograms`
        for use later (plotting, outputting)

        :return: None
        """

        # convert to overall distribution of all box densities and find "middle density"
        all_box_densities = []
        for box in self.total_densities.keys():
            all_box_densities += list(self.total_densities[box].values())
        all_box_densities = np.array(all_box_densities)
        middle_density = all_box_densities.mean()

        # partition all densities into low and high regions
        density_vals = {key: {} for key in self.region_names}
        for box in self.total_densities.keys():
            for cycle, value in self.total_densities[box].items():
                if value < middle_density:
                    region = self.region_names[0]
                else:
                    region = self.region_names[1]
                key = '{}-{}'.format(cycle, box)  # keep box information in data for later
                density_vals[region][key] = value

        histogram_data = {key: {} for key in self.region_names}
        # iterate through each region and do fitting
        for region in self.region_names:
            values_for_binning = density_vals[region]

            vals = np.array(list(values_for_binning.values()))
            counts, edges = np.histogram(vals,
                                         bins=self.number_of_bins,
                                         density=False)
            total_count = sum(counts)
            # normalize so that sum of all probabilities is 1
            probabilities = [i / total_count for i in counts]

            # nodes are in between edges
            nodes = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
            histogram_data[region]['nodes'] = nodes[:]
            histogram_data[region]['vals'] = probabilities[:]

            # find subset of data that is above threshold
            indices, = np.where(probabilities >= self.fraction_peak_maximum_density * np.max(probabilities))
            x = [nodes[i] for i in indices]
            y = [probabilities[i] for i in indices]
            histogram_data[region]['nodes above threshold'] = x[:]
            histogram_data[region]['vals above threshold'] = y[:]

            (params, covariance, min_gauss, max_gauss) = self.fit_gaussian(x, y)

            # find densities, pressures, compressibilities in each region
            for key in values_for_binning.keys():
                cycle, box = key.split('-')
                cycle = int(cycle)
                self.region_data['density'][region][key] = values_for_binning[key]
                if cycle in self.pressures[box].keys():
                    self.region_data['pressure'][region][key] = self.pressures[box][cycle]
                if cycle in self.compressibilities[box].keys():
                    self.region_data['compressibility'][region][key] = self.compressibilities[box][cycle]
                if max_gauss >= values_for_binning[key] >= min_gauss:
                    self.region_data_above_threshold['density'][region][key] = values_for_binning[key]
                    if cycle in self.pressures[box].keys():
                        self.region_data_above_threshold['pressure'][region][key] = self.pressures[box][cycle]
                    if cycle in self.compressibilities[box].keys():
                        self.region_data_above_threshold['compressibility'][region][key] = self.compressibilities[box][cycle]

            histogram_data[region]['raw data'] = list(self.region_data_above_threshold['density'][region].values())

            # save gaussian fit results
            self.gaussian_fit_results[region] = {
                'parameters': params,
                'covariance': covariance,
                'min_gauss': min_gauss,
                'max_gauss': max_gauss
            }
        self.histograms['density'] = histogram_data

    def histogram_P_Z_data(self, data):
        """Obtain histogram(s) of pressure or compressibility data within region(s) or box(es)

        :param data:
        :type data: dict
        :return: results of histograms
        :rtype: dict
        """
        hist_data = {}
        for key in data.keys():
            hist_data[key] = {}
            data_vals = np.array(list(data[key].values()))
            hist, edges = np.histogram(data_vals,
                                       bins=self.number_of_bins,
                                       density=False)
            bin_means = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
            total_count = sum(hist)
            probabilities = [i/total_count for i in hist]
            hist_data[key] = {
                'nodes': bin_means,
                'vals': probabilities,
                'raw data': list(data[key].values())
            }
        return hist_data

    def calculate_P_Z(self):
        """Calculate pressure and compressibilities for two different types of distributions

        Store data of mean and standard deviation of raw data values in
        :py:attr:`~histograms.density_histogram.DPD.histograms`
        within
        :py:attr:`~histograms.density_histogram.DPD.data`

        :return: None
        """
        self.histograms['pressure above threshold'] = self.histogram_P_Z_data(self.region_data_above_threshold['pressure'])
        self.histograms['compressibility above threshold'] = self.histogram_P_Z_data(self.region_data_above_threshold['compressibility'])
        self.histograms['pressure'] = self.histogram_P_Z_data(self.region_data['pressure'])
        self.histograms['compressibility'] = self.histogram_P_Z_data(self.region_data['compressibility'])
        self.store_data()

    def store_data(self):
        """
        Store data of mean and standard deviation of raw data values in
        :py:attr:`~histograms.density_histogram.DPD.histograms`
        within
        :py:attr:`~histograms.density_histogram.DPD.data`

        :return: None
        """
        self.data = {}
        for data_type in self.histograms.keys():
            for peak, vals in self.histograms[data_type].items():
                if peak not in self.data.keys():
                    self.data[peak] = {}
                mean = np.mean(vals['raw data'])
                stdev = np.std(vals['raw data'])
                self.data[peak][data_type] = {
                    'mean': mean,
                    'stdev': stdev
                }

    def output_results(self):
        """
        Output information stored in
        :py:attr:`~histograms.density_histogram.DPD.data`

        Change to g/mL units if was desired during instantiation

        Output results either by printing or saving to csv file,
            as requested by input arguments

        :return: None
        """
        if (self.unit == "g/mL"):
            self.units["density"] = "g/mL"
            low_key, high_key = self.region_names
            self.data[high_key]['density']['mean'] *= (1 / N_av) * (1 / 1e-21) * (self.molecular_weight)
            self.data[high_key]['density']['stdev'] *= (1 / N_av) * (1 / 1e-21) * (self.molecular_weight)
            self.data[low_key]['density']['mean'] *= (1 / N_av) * (1 / 1e-21) * (self.molecular_weight)
            self.data[low_key]['density']['stdev'] *= (1 / N_av) * (1 / 1e-21) * (self.molecular_weight)
        if not self.output_file_name:
            import pprint
            pprint.pprint(self.data)
        else:
            with open(self.output_path + self.output_file_name + '_histogram_averages.csv', 'w') as f:
                f.write('Property,Units,High Density Mean,High Density Standard Deviation,'
                        'Low Density Mean,Low Density Standard Deviation\n')
                low_key, high_key = self.region_names
                for property in 'pressure above threshold', 'compressibility above threshold', 'density':
                    f.write(
                        '%s,%s,%e,%e,%e,%e\n' % (
                            property, self.units[property.split()[0]],
                            self.data[high_key][property]['mean'], self.data[high_key][property]['stdev'],
                            self.data[low_key][property]['mean'], self.data[low_key][property]['stdev']
                        )
                    )
        self.output_blocks()
        self.output_histograms()
        self.output_gaussian_fit()
        self.save_region_specific_trajectory()

    def output_histograms(self):
        with open(self.output_path+self.output_file_name+'_histograms.json', 'w') as f:
            json.dump(self.histograms, f)

    def output_gaussian_fit(self):
        with open(self.output_path+self.output_file_name+'_gaussian_results.json', 'w') as f:
            json.dump(self.gaussian_fit_results, f)

    def save_region_specific_trajectory(self):
        with open(self.output_path+self.output_file_name+'_region_data_above_threshold.json', 'w') as f:
            json.dump(self.region_data_above_threshold, f)
        with open(self.output_path+self.output_file_name+'_region_data.json', 'w') as f:
            json.dump(self.region_data, f)

    def decide_if_two_plots(self, data):
        means = (data[self.region_names[0]]['parameters'][-2],
                 data[self.region_names[1]]['parameters'][-2])
        if max(means) > 10. * min(means):
            return True
        else:
            return False

    def plot(self, data_type):
        self.figure_number += 1
        fig = plt.figure(self.figure_number)
        axes = [fig.add_subplot(111)]
        two_plots = False
        if data_type == 'density':
            two_plots = self.decide_if_two_plots(self.gaussian_fit_results)
        if two_plots:
            self.figure_number += 1
            fig2 = plt.figure(self.figure_number)
            axes.append(fig2.add_subplot(111))
        colors = {
            self.region_names[0]: {'nodes': 'C0', 'fit': 'C1', 'mean': 'C4'},
            self.region_names[1]: {'nodes': 'C2', 'fit': 'C3', 'mean': 'C5'},
            'box 1': {'nodes': 'C6', 'mean': 'C7'},
            'box 2': {'nodes': 'C8', 'mean': 'C9'}
        }
        for i, key in enumerate(self.histograms[data_type].keys()):
            if i < len(axes):
                ax = axes[i]
            else:
                ax = axes[-1]
            values = self.histograms[data_type][key]
            nodes = values['nodes']
            vals = values['vals']
            ax.plot(nodes, vals, '.', color=colors[key]['nodes'], label='nodes, %s' % key)
            if 'parameters' in values.keys():
                # density
                gauss_x = np.linspace(values['min_gauss'], values['max_gauss'])
                gauss_y = gaussian(gauss_x, *values['parameters'])
                ax.plot(gauss_x, gauss_y, ls='solid', color=colors[key]['fit'], label='fit, %s' % key)
            else:
                # pressure, compressibility
                mean_value = np.mean(values['raw data'])
                ax.plot(mean_value, np.max(vals), 'x', color=colors[key]['nodes'],
                        label='mean, %s' % key)
        if data_type + ' all' in self.histograms.keys():
            prob_distribution_and_fit = self.histograms[data_type + ' all']
            for i, key in enumerate(prob_distribution_and_fit.keys()):
                values = prob_distribution_and_fit[key]
                nodes = values['nodes']
                vals = values['vals']
                if i < len(axes):
                    ax = axes[i]
                else:
                    ax = axes[-1]
                ax.plot(nodes, vals, '.', color=colors[key]['nodes'], label='nodes, %s' % key)
                mean_value = np.mean(values['raw data'])
                ax.plot(mean_value, np.max(vals), 's', color=colors[key]['nodes'],
                        label='mean, %s' % key)
        for ax in axes:
            ax.set_ylabel('Probability')
            ax.set_xlabel(data_type)
            ax.legend()
        if self.output_file_name:
            if two_plots:
                fig.savefig(self.output_path + self.output_file_name + '_1_' + data_type + '.png')
                fig2.savefig(self.output_path + self.output_file_name + '_2_' + data_type + '.png')
            else:
                fig.savefig(self.output_path + self.output_file_name + '_' + data_type + '.png')
        else:
            plt.show()
        plt.close(fig)
        if two_plots:
            plt.close(fig2)

    def main(self):
        """
        Main method of :class:
        First, calls :py:func:`~histogram.density_histogram.DPD.choose_output_file_path`
        Then, calls :py:func:`~histogram.density_histogram.DPD.initialize_data`
        The peaks are fit using
        :py:func:`~histogram.density_histogram.DPD.fit_density_peaks`
        and the corresponding pressure and compressibilities are obtained from
        :py:func:`~histogram.density_histogram.DPD.calculate_P_Z`
        and then plot results if :py:attr:`~histogram.density_histogram.DPD.make_plot` == 'Yes'
        and then output results with
        :py:func:`~hisogram.density_histogram.DPD.output_results`
        :return: None
        """
        self.choose_output_file_path()
        self.initialize_data()
        self.fit_density_peaks()
        self.calculate_P_Z()
        if self.make_plot == 'Yes':
            for key in 'pressure', 'compressibility', 'density':
                self.plot(key)
        self.output_results()


def parse_arguments():
    """Command-line parser for input data

    :return: arguments to pass to either :py:class:`Blocks` or :py:class:`DPD` during instantiation
    :rtype: dict
    """
    import argparse
    parser = argparse.ArgumentParser(description='Read trajectory file(s). Calculate densities'
                                                 ' of two different boxes using the histogram (gaussian) approach.'
                                                 ' Calculate vapor pressure (and compressibility factor) using data'
                                                 ' for the high-density box and low-density box')
    parser.add_argument('-fi', '--files', help='list of trajectory files for ONE independent simulation',
                        type=str, nargs='+')
    parser.add_argument('-fd', '--fraction_peak_maximum_density',
                        help='Fraction of peak maximum to consider in gaussian fitting of densities'
                             '(defaults to 0.75; see J. Chem. Phys. 143, 114113 (2015) for recommendations)',
                        type=float, default=0.75)
    parser.add_argument('-n', '--number_of_bins', help='Number of bins within 2 stdevs of each density peak'
                                                       '(defaults to 100)',
                        type=int, default=100)
    parser.add_argument('-T', '--temperature', help='Temperature of simulation [ Kelvin units ]', type=float)
    parser.add_argument('-o', '--output_file_name',
                        help='Name of output file (optional, used for results file and plot)', type=str,
                        default='output')
    parser.add_argument('-p', '--plot', help='whether to plot results', type=str,
                        choices=['Yes', 'No'], default='No')
    parser.add_argument('-b', '--blocks_only', help='whether or not to only do block averages only',
                        choices=['Yes', 'No'], default='No')
    parser.add_argument('-u', '--d_units', help='Choice for density units', type=str,
                        choices=['chains/nm**3', 'g/mL'], default='chains/nm**3')
    args = vars(parser.parse_args())
    return args
