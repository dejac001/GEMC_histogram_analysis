from __future__ import division, print_function


def read_fort12(file_name, cycle_start):
    """Read fort.12 (trajectory) file for MCCCS-MN code output file style.
    Most inner data type is {cycle: value}. This is helpful for finding out
    what cycle a certain value occured at, etc.

    :param file_name:   name of file to open and read
    :type file_name: str
    :param cycle_start: cycle number to start at (for labeling data explicitly at each cycle)
    :type cycle_start: int
    :returns:  (molecular_weights, number_of_molecules, box_volume, box_pressure, iratp, cycle_counter)
    """
    f = open(file_name)
    line_number = 1
    first_line = next(f)
    (number_of_steps_planned, iratp, number_of_boxes, number_of_moltypes), MWs = map(int, first_line.split()[
                                                                                           :4]), first_line.split()[4:]
    assert len(MWs) == number_of_moltypes, 'Inconsistent number of molecule types (%i) and weights (%i)' % (
        number_of_moltypes, len(MWs))
    cycle_counter = cycle_start
    number_of_molecules = {
        'molecule %i' % i: {'box %i' % j: {} for j in range(1, number_of_boxes + 1)}
        for i in range(1, number_of_moltypes + 1)
    }
    box_pressure = {'box %i' % i: {} for i in range(1, number_of_boxes + 1)}
    box_internal_energy = {'box %i' % i: {} for i in range(1, number_of_boxes + 1)}
    box_volume = {'box %i' % i: {} for i in range(1, number_of_boxes + 1)}
    molecular_weights = {'molecule %i' % i: float(MWs[i - 1]) for i in range(1, number_of_moltypes + 1)}

    for line in f:  # starting on 2nd line
        line_number += 1
        box = line_number % number_of_boxes + 1
        box_name = 'box %i' % box
        if box % number_of_boxes == 1:
            cycle_counter += 1
        my_split = line.split()
        offset = len(my_split) - number_of_moltypes - 1
        for mol in range(1, number_of_moltypes + 1):
            number_of_molecules['molecule %i' % mol][box_name][cycle_counter] = int(my_split[offset + mol])
        # TODO: add implementation for nonorthorhombic boxes
        volume = float(my_split[0]) * float(my_split[1]) * float(my_split[2]) / 1000.  # convert to nm**3
        box_volume[box_name][cycle_counter] = volume
        box_internal_energy[box_name][cycle_counter] = float(my_split[3])
        if cycle_counter > cycle_start and (cycle_counter-cycle_start) % iratp == 0:
            #  pressure was calculated on this run
            box_pressure[box_name][cycle_counter] = float(my_split[4])
    return (
        molecular_weights,
        number_of_molecules,
        box_volume,
        box_pressure,
        iratp,
        cycle_counter
    )
