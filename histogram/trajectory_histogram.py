import matplotlib.pyplot as plt

from histogram import reader
from histogram.density_histogram import Blocks, update_data_box, update_data_molty


class Trajectory(Blocks):
    def __init__(self, **kwargs):
        Blocks.__init__(self, **kwargs)
        self.num_molecules = None
        self.number_of_bins = kwargs['number_of_bins']
        self.box_volumes = None
        self.total_volume = None
        self.total_number_of_molecules = None

    def initialize_data(self):
        all_data = {key: {} for key in ['N', 'V']}
        for file in self.files:
            f_nmolty, f_volume, f_pressure, f_N_with_p, f_V_with_p, iratp = reader.read_fort12(file)
            if not self.iratp:
                self.iratp = iratp
            update_data_molty(all_data['N'], f_nmolty)
            update_data_box(all_data['V'], f_volume)
        self.box_volumes = all_data['V']
        self.boxes = list(i for i in self.box_volumes.keys())
        self.total_number_of_molecules = sum(
            all_data['N'][i][box][0] for i in all_data['N'].keys() for box in all_data['V'].keys()
        )
        self.num_molecules = {box: [] for box in self.box_volumes}
        for box in self.boxes:
            for i in range(len(all_data['N']['molecule 1'][box])):
                self.num_molecules[box].append(
                    sum(all_data['N'][mol][box][i] for mol in all_data['N'].keys())
                )
        # TODO: check that total volume does not change
        self.total_volume = sum(
            all_data['V'][box][0] for box in all_data['V'].keys()
        )  # nm**3
        self.box_volumes = {box:
                                [i / self.total_volume for i in self.box_volumes[box]]
                            for box in all_data['V'].keys()}
        self.num_molecules = {box:
                                  [i / self.total_number_of_molecules for i in self.num_molecules[box]]
                              for box in all_data['V'].keys()}

    def plot_probability_histogram(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        Ni_over_Nt = []
        Vi_over_Vt = []
        for box in self.boxes:
            Ni_over_Nt = Ni_over_Nt + self.num_molecules[box]
            Vi_over_Vt = Vi_over_Vt + self.box_volumes[box]
        hist, x_edges, y_edges, im = ax.hist2d(
            Vi_over_Vt, Ni_over_Nt, bins=self.number_of_bins, normed=True,
            range=[[0., 1.], [0., 1.]],
            cmin=1,  # white space if count is less than 1
        )
        cax = fig.add_axes([0.77, 0.11, 0.05, 0.87])
        cbar = fig.colorbar(im, cax=cax)
        from plotting_util import set_x_ticks, set_y_ticks
        set_x_ticks(ax, [0, 0.2, 0.4, 0.6, 0.8, 1.0])
        set_y_ticks(ax, [0, 0.2, 0.4, 0.6, 0.8, 1.0])
        cbar.set_label('Probability Distribution')
        cax.tick_params(which='both', direction='out')
        ax.set_ylabel('$N_i/N_{\mathrm{total}}$')
        ax.set_xlabel('$V_i/V_{\mathrm{total}}$')
        plt.subplots_adjust(
            left=0.13, right=0.72
        )
        if self.output_file_name:
            fig.savefig(self.output_path + self.output_file_name + '_' +
                        'sampling-trajectory.png')
        else:
            plt.show()
        plt.close(fig)

    def main(self):
        """Main workflow for trajectory

        :return: None
        """
        self.choose_output_file_path()
        self.initialize_data()
        self.plot_probability_histogram()
