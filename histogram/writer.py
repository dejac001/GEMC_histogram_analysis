def histogram_bins_and_values(output_prefix, data):
    for key in data.keys():
        factor = 1.
        if 'pressure' in key:
            factor = 1./1000
        for region, val in data[key].items():
            y_scale = 1.
            if 'above threshold' in key:
                y_scale = sum(data['density'][region]['vals above threshold'])/sum(data['density'][region]['vals'])
            with open(output_prefix + '%s_%s_hist.dat' % (key, region.replace(' ', '-')), 'w') as f:
                for i, j in zip(val['nodes'], val['vals']):
                    f.write('%e %e\n' % (i*factor, j*y_scale))
                if key == 'density':
                    # special density histogram above peak
                    with open(output_prefix + 'density_peak_%s_hist.dat' % region.replace(' ', '-'), 'w') as f:
                        for i, j in zip(val['nodes above threshold'],
                                        val['vals above threshold']):
                            f.write('%e %e\n' % (i, j))


def gaussian_function(output_prefix, data):
    import numpy as np
    from histogram.density_histogram import gaussian
    for region, val in data.items():
        x = np.linspace(val['min_gauss'], val['max_gauss'], 100)
        y = gaussian(x, *val['parameters'])
        with open(output_prefix + 'gaussian_%s.dat' % region, 'w') as f:
            for i, j in zip(x, y):
                f.write('%e %e\n' % (i, j))

