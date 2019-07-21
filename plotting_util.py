def set_y_ticks(my_ax, ticks):
    my_ax.set_yticks(ticks)
    y_minor = [(ticks[i] + ticks[i + 1]) / 2. for i in range(len(ticks) - 1)]
    my_ax.set_yticks(y_minor, minor=True)
    my_ax.set_ylim([min(ticks), max(ticks)])


def set_x_ticks(my_ax, ticks):
    my_ax.set_xticks(ticks)
    x_minor = [(ticks[i] + ticks[i + 1]) / 2. for i in range(len(ticks) - 1)]
    my_ax.set_xticks(x_minor, minor=True)
    my_ax.set_xlim([min(ticks), max(ticks)])
