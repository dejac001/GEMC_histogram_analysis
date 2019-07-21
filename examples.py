"""
example workflows for using the histogram code. To run the examples,
call the functions in this file
"""


from analysis_factory import RunAllIndepSims, average


def butane(T):
    """
    analyze the data for butane at T K, located in the example_data/butane_XXXK directory
        where XXX is the temperature
    """
    results = RunAllIndepSims(temperature=T, n_bins=100)
    indep_sim_list = list(map(str, range(1, 9)))
    for sim in indep_sim_list:
        results.run_one_indep_sim('example_data/butane_%iK/%s/' % (T, sim), plot='Yes')
    average('example_data/butane_%iK' % T, indep_sim_list)


def butane_400K():
    return butane(400)


def butane_264K():
    return butane(264)


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


def pentanal_360K_short(**kwargs):
    indep_sims = list('par'+i for i in  map(str, range(1, 9)))
    pentanal(indep_sims, T=360, **kwargs)


def pentanal_long(temperature, number_of_bins=50, average_fmt='-long'):
    indep_sims = ['par%s' % i for i in map(str, range(1, 17))]
    pentanal(indep_sims, **{'T': temperature, 'number_of_bins': number_of_bins,
                            'average_fmt': '-long'})


def example_bin_width():
    """
    analyze effect of number of bins chosen on
        the data for pentanal at 360K, located in the example_data/pentanal_360K directory
    """
    for n_bins in [20, 40, 100]:
        pentanal_360K_short(number_of_bins=n_bins, average_fmt='-%i-bins' % n_bins)

