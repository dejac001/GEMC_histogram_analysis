"""
example workflows for using the histogram code. To run the examples,
call the functions in this file
"""


from analysis_factory import RunAllIndepSims, average


def pentanal(list_of_indep_sims, T=0, number_of_bins=0, average_fmt=''):
    """
    analyze the data for pentanal at a given temperature=T, located in the example_data/pentanal_XXXK directory
        where XXX is the temperature

    :param list_of_indep_sims: independent simulations for which to average over
    :param T: temperature of simulation in K
    :param number_of_bins: number of bins to analyze over for each region
    :param average_fmt: subscript for identifier in data output
    :type average_fmt: str
    :type number_of_bins: int
    :type T: float
    :type list_of_indep_sims: list
    """
    results = RunAllIndepSims(temperature=T, n_bins=number_of_bins)
    for dir in list_of_indep_sims:
        results.run_one_indep_sim('example_data/pentanal_%iK/%s' % (T, dir), plot='Yes',
                                  output_file_name='output' + average_fmt)
    average('example_data/pentanal_%iK' % T, list_of_indep_sims, fmt=average_fmt)


def pentanal_long(temperature, number_of_bins=50, average_fmt='-long'):
    indep_sims = ['par%s' % i for i in map(str, range(1, 17))]
    pentanal(indep_sims, **{'T': temperature, 'number_of_bins': number_of_bins,
                            'average_fmt': '-long'})

