from histogram.density_histogram import parse_arguments, Blocks, DPD


def run():
    """
    Perform histogram analysis. Main function of run_histogram.py.
    :return: None
    """
    args = parse_arguments()
    if args['blocks_only'] == 'Yes':
        instance = Blocks(**args)
    else:
        instance = DPD(**args)
    instance.main()


if __name__ == '__main__':
    run()
