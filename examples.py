"""
example workflows for using the histogram code. To run the examples,
call the functions in this file
"""


from analysis_factory import RunAllIndepSims, average


def pentanal_510K(number_of_bins=0):
    """
    analyze the data for pentanal for the first independent simulation at 510 K

    :param number_of_bins: number of bins to analyze over for each region
    :type number_of_bins: int
    """
    T = 510
    results = RunAllIndepSims(temperature=T, n_bins=number_of_bins)
    dir = 'par1'
    results.run_one_indep_sim('example_data/pentanal_%iK/%s' % (T, dir), plot='Yes',
                                  output_file_name='output')
    average('example_data/pentanal_%iK' % T, ['par1'])

