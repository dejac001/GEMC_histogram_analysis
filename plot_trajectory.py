from histogram.density_histogram import parse_arguments
from histogram.trajectory_histogram import Trajectory


def example_butane():
    args = parse_arguments()
    args['files'] = [
        'example_data/butane_400K/1/fort12.prod-%i' % i for i in range(1, 4)
    ]
    args['temperature'] = 400.
    args['plot'] = 'Yes'
    instance = Trajectory(**args)
    instance.main()
    return instance


def example_pentanal():
    args = parse_arguments()
    args['files'] = [
        'example_data/pentanal_510K/par1/fort12.prod%i' % i for i in range(1, 20)
    ]
    args['temperature'] = 510.
    args['plot'] = 'Yes'
    instance = Trajectory(**args)
    instance.main()
    return instance


def run_examples():
    example_pentanal()
    example_butane()


if __name__ == '__main__':
    # args = parse_arguments()
    # instance = Trajectory(**args)
    # instance.main()
    run_examples()
