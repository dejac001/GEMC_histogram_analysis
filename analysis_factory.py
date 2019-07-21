import os

from histogram.density_histogram import DPD


class RunAllIndepSims:
    def __init__(self, temperature, n_bins=100):
        self.T = temperature
        self.instances = []
        self.number_of_figures_plotted = 0
        self.n_bins = n_bins

    def run_one_indep_sim(self, dir_prefix, plot='Yes', output_file_name='output',
                          fraction_peak_maximum_density=0.75):
        args = {'fraction_peak_maximum_density': fraction_peak_maximum_density, 'temperature': self.T, 'plot': plot,
                'number_of_bins': self.n_bins, 'output_file_name': output_file_name,
                'd_units': 'default'}
        fort_12_files = []
        for file in os.listdir(dir_prefix):
            # try to automatically find all fort.12 files; might not work for everyone
            # might not work for everyone
            if 'fort12' in file or ('fort' in file and '12' in file):
                fort_12_files.append(file)
        args['files'] = [dir_prefix + '/' + i for i in fort_12_files]
        instance = DPD(**args)
        instance.figure_number = self.number_of_figures_plotted  # update global number of figs plotted
        instance.main()
        self.add_indep_sim(instance)
        self.number_of_figures_plotted += instance.figure_number
        self.average_all_indep_sims()

    def add_indep_sim(self, instance):
        # TODO: add new data from indep sim to class
        #   dont add all data, then the amount of memory used can get excessive
        pass

    def average_all_indep_sims(self):
        # TODO: implement this by iterating over some data type created in self.add_indep sim
        pass


def csv_to_dict(file_name):
    import csv
    with open(file_name) as f:
        data = csv.reader(f)
        ordered_keys = next(data)
        my_dict = {key: [] for key in ordered_keys}
        for line_data in data:
            assert len(line_data) == len(ordered_keys), 'Inconsistent column formatting'
            for key, val in zip(ordered_keys, line_data):
                my_dict[key].append(val)
    return my_dict


def average(dir, indep_dirs, fmt=''):
    if dir[-1] != '/':
        dir = dir + '/'
    for seed in indep_dirs:
        hist_file = dir + '/%s/output%s_histogram_averages.csv' % (seed, fmt)
        block_file = dir + '/%s/output%s_block_averages.csv' % (seed, fmt)
        hist_data = csv_to_dict(hist_file)
        block_data = csv_to_dict(block_file)
        if seed == indep_dirs[0]:
            averages = {prop: {} for prop in hist_data['Property']}
        for key in hist_data.keys():
            if 'Mean' in key:
                for index, property in enumerate(hist_data['Property']):
                    if seed == indep_dirs[0]:
                        averages[property][key] = []
                    averages[property][key].append(float(hist_data[key][index]))
        for key in block_data.keys():
            if 'Mean' in key:
                for index, property in enumerate(hist_data['Property']):
                    if seed == indep_dirs[0]:
                        averages[property][key] = []
                    averages[property][key].append(float(block_data[key][index]))

    import numpy as np
    with open(dir + 'results%s.csv' % fmt, 'w') as f:
        f.write('%-16s,%-21s,%-12s,%-22s\n' %
                ('Property', 'Key', 'Mean', 'Standard Error of Mean'))
        for property in averages.keys():
            for key, val in averages[property].items():
                f.write('%-16s,%-21s,%-12f,%-22f\n' %
                        (property, key, np.mean(np.array(val)), np.std(np.array(val), ddof=1) / np.sqrt(len(val)))
                        )
